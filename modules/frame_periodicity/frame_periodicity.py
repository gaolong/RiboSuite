#!/usr/bin/env python3

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../../"))
sys.path.insert(0, ROOT_DIR)

import argparse
from collections import defaultdict

import pysam
import pandas as pd

# --- plotting (HPC / headless safe) ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --------------------------------------------------
# shared utilities
# --------------------------------------------------
from bin.gtf_models import load_tx_models
from bin.psite_utils import (
    load_offsets,
    get_psite,
    in_interval,
    genomic_to_tx_spliced,
)

# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Frame periodicity QC for 5'UTR, CDS and 3'UTR (CDS-anchored frames)"
    )
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--psite_offsets", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min_cds_len", type=int, default=0)
    ap.add_argument("--top_tx", type=int, default=15000)
    ap.add_argument("--mapq", type=int, default=100)

    # Keep for backward-compat; not used with CDS-anchored models
    ap.add_argument(
        "--use_phase",
        action="store_true",
        help="(ignored) Phase is not used in CDS-anchored model periodicity."
    )
    return ap.parse_args()


# --------------------------------------------------
# choose the longest tx for each gene
# --------------------------------------------------
def keep_longest_tx_per_gene(models):
    """
    Keep one transcript per gene: the one with the longest CDS length.
    """
    best = {}
    for m in models:
        gene = m["gene"]
        cds_len = m.get("cds_len", 0)
        if gene not in best or cds_len > best[gene]["cds_len"]:
            best[gene] = {"model": m, "cds_len": cds_len}
    return [v["model"] for v in best.values()]


# --------------------------------------------------
# region assignment in tx coords (tx coords increase with genomic coordinate)
# --------------------------------------------------
def assign_region_tx(tx_pos, strand, cds_lo_tx, cds_hi_tx, utr_tx):
    """
    Returns region in {"5UTR","CDS","3UTR"} or None if outside modeled exons.
    Uses utr_tx intervals (already exon - CDS).
    """
    in_cds = (cds_lo_tx <= tx_pos < cds_hi_tx)
    if in_cds:
        return "CDS"

    if not in_interval(tx_pos, utr_tx):
        return None

    # Plus-strand: 5' is lower tx coords
    if strand == "+":
        return "5UTR" if tx_pos < cds_lo_tx else "3UTR"

    # Minus-strand: 5' is higher tx coords (because tx coords follow genomic increasing)
    return "5UTR" if tx_pos >= cds_hi_tx else "3UTR"


# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()

    offsets = load_offsets(args.psite_offsets)
    if not offsets:
        raise RuntimeError("No offsets loaded")
    max_offset = max(offsets.values())

    bam = pysam.AlignmentFile(args.bam, "rb")

    # Load models from shared bin/gtf_models.py
    models = load_tx_models(
        args.gtf,
        min_cds_len=args.min_cds_len,
        protein_coding_only=True,
    )
    models = keep_longest_tx_per_gene(models)

    if not models:
        raise RuntimeError("No valid transcripts found")

    # --------------------------------------------------
    # PASS 1: select top transcripts by CDS P-site counts
    # --------------------------------------------------
    tx_counts = defaultdict(int)

    for m in models:
        chrom, strand = m["chrom"], m["strand"]
        tx_id = (m["gene"], m["tx"])

        fs = max(0, m["fetch_start"] - max_offset)
        fe = m["fetch_end"] + max_offset

        for read in bam.fetch(chrom, fs, fe):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq:
                continue

            psite, rl = get_psite(read, offsets, strand=strand)
            if psite is None:
                continue

            tx_pos = genomic_to_tx_spliced(psite, m["exons"])
            if tx_pos is None:
                continue

            if m["cds_lo_tx"] <= tx_pos < m["cds_hi_tx"]:
                tx_counts[tx_id] += 1

    top_tx = {
        tx for tx, _ in sorted(tx_counts.items(), key=lambda x: x[1], reverse=True)[:args.top_tx]
    }

    # --------------------------------------------------
    # PASS 2: frame periodicity
    # --------------------------------------------------
    counts = defaultdict(int)

    for m in models:
        tx_id = (m["gene"], m["tx"])
        if tx_id not in top_tx:
            continue

        chrom, strand = m["chrom"], m["strand"]
        fs = max(0, m["fetch_start"] - max_offset)
        fe = m["fetch_end"] + max_offset

        anchor_tx = m["anchor_tx"]
        cds_lo_tx = m["cds_lo_tx"]
        cds_hi_tx = m["cds_hi_tx"]
        utr_tx = m["utr_tx"]

        for read in bam.fetch(chrom, fs, fe):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq:
                continue

            psite, rl = get_psite(read, offsets, strand=strand)
            if psite is None:
                continue

            tx_pos = genomic_to_tx_spliced(psite, m["exons"])
            if tx_pos is None:
                continue

            region = assign_region_tx(tx_pos, strand, cds_lo_tx, cds_hi_tx, utr_tx)
            if region is None:
                continue

            # Frame (CDS anchored)
            if strand == "+":
                frame = (tx_pos - anchor_tx) % 3
            else:
                frame = (anchor_tx - tx_pos) % 3

            counts[(region, rl, frame)] += 1

    # --------------------------------------------------
    # output
    # --------------------------------------------------
    df = pd.DataFrame([
        dict(sample=args.sample, region=r, read_length=rl, frame=f, count=c)
        for (r, rl, f), c in counts.items()
    ])

    out_tsv = f"{args.sample}.periodicity.by_region_by_length.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[periodicity] written {out_tsv}")

    if df.empty:
        return

    df["fraction"] = df["count"] / df.groupby(["region", "read_length"])["count"].transform("sum")

    # --------------------------------------------------
    # plotting
    # --------------------------------------------------
    regions = ["5UTR", "CDS", "3UTR"]
    titles = ["5′ UTR", "CDS", "3′ UTR"]

    cmap = LinearSegmentedColormap.from_list("white_red", ["#ffffff", "#ff0000"])
    lengths = sorted(df["read_length"].unique())

    fig, axes = plt.subplots(
        ncols=3,
        figsize=(11, max(4, len(lengths) * 0.35)),
        constrained_layout=True
    )

    for ax, region, title in zip(axes, regions, titles):
        sub = df[df["region"] == region]

        mat = (
            sub.pivot_table(
                index="read_length",
                columns="frame",
                values="fraction",
                fill_value=0.0
            )
            .reindex(lengths, fill_value=0.0)
        )

        im = ax.imshow(
            mat.values,
            aspect="auto",
            cmap=cmap,
            vmin=0,
            vmax=1,
            origin="upper"
        )

        ax.set_title(title)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["Frame 0", "Frame 1", "Frame 2"])

        ax.set_yticks(range(len(lengths)))
        if ax is axes[0]:
            ax.set_yticklabels(lengths)
            ax.set_ylabel("Read length")
        else:
            ax.set_yticklabels([])

    fig.colorbar(im, ax=axes, fraction=0.035, pad=0.02).set_label("Fraction")
    fig.suptitle(f"{args.sample} frame periodicity")

    out_png = f"{args.sample}.periodicity.by_region_by_length.heatmap.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"[periodicity] written {out_png}")


if __name__ == "__main__":
    main()