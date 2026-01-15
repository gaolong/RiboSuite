#!/usr/bin/env python3

import pysam
import argparse
import pandas as pd
from collections import defaultdict
import re

# --- plotting (HPC / headless safe) ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, FixedFormatter

# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Frame periodicity QC for 5'UTR, CDS and 3'UTR "
                    "(ribowaltz-style, CDS-anchored frames)"
    )
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min_cds_len", type=int, default=150)
    ap.add_argument("--top_tx", type=int, default=1000)
    ap.add_argument("--mapq", type=int, default=10)
    return ap.parse_args()


# --------------------------------------------------
# offsets
# --------------------------------------------------
def load_offsets(path):
    df = pd.read_csv(path, sep="\t")
    col = "psite_offset" if "psite_offset" in df.columns else "offset"
    return dict(zip(df["read_length"].astype(int), df[col].astype(int)))


# --------------------------------------------------
# GTF helpers
# --------------------------------------------------
def parse_attrs(attr):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr))


def load_start_codons(gtf):
    start_tx = set()
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, *_ , attr = line.rstrip().split("\t")
            if feature != "start_codon":
                continue
            attrs = parse_attrs(attr)
            if "gene_id" in attrs and "transcript_id" in attrs:
                start_tx.add((attrs["gene_id"], attrs["transcript_id"]))
    return start_tx


# --------------------------------------------------
# genomic → spliced transcript coordinate
# --------------------------------------------------
def genomic_to_tx_spliced(pos, exons, strand):
    if strand == "+":
        for s, e, tx_s in exons:
            if s <= pos < e:
                return tx_s + (pos - s)
    else:
        for s, e, tx_s in reversed(exons):
            if s <= pos < e:
                return tx_s + (e - pos - 1)
    return None


# --------------------------------------------------
# load transcript models (exon-based)
# --------------------------------------------------
def load_tx_models(gtf, min_cds_len):
    tx_exons = defaultdict(list)
    tx_cds = defaultdict(list)
    gene_type = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            chrom, _, feature, start, end, _, strand, phase, attr = line.rstrip().split("\t")
            attrs = parse_attrs(attr)

            if feature == "gene":
                gene_type[attrs["gene_id"]] = (
                    attrs.get("gene_type") or attrs.get("gene_biotype")
                )
                continue

            gene = attrs.get("gene_id")
            tx = attrs.get("transcript_id")
            if not gene or not tx:
                continue
            if gene_type.get(gene) != "protein_coding":
                continue

            start0 = int(start) - 1
            end0 = int(end)

            if feature == "exon":
                tx_exons[(gene, tx, chrom, strand)].append((start0, end0))
            elif feature == "CDS":
                tx_cds[(gene, tx, chrom, strand)].append((start0, end0))

    models = []

    for key, exons in tx_exons.items():
        gene, tx, chrom, strand = key
        if key not in tx_cds:
            continue

        exons.sort(key=lambda x: x[0])
        cds_intervals = tx_cds[key]
        cds_len = sum(e - s for s, e in cds_intervals)
        if cds_len < min_cds_len:
            continue

        # build spliced transcript coordinates
        tx_pos = []
        cursor = 0
        for s, e in exons:
            tx_pos.append((s, e, cursor))
            cursor += (e - s)
        tx_len = cursor

        # CDS boundaries in genomic coords
        cds_start_g = min(s for s, e in cds_intervals)
        cds_end_g = max(e for s, e in cds_intervals)

        cds_start_tx = genomic_to_tx_spliced(cds_start_g, tx_pos, strand)
        cds_end_tx = genomic_to_tx_spliced(cds_end_g - 1, tx_pos, strand) + 1

        models.append({
            "gene": gene,
            "tx": tx,
            "chrom": chrom,
            "strand": strand,
            "exons": tx_pos,
            "tx_len": tx_len,
            "cds_start_tx": cds_start_tx,
            "cds_end_tx": cds_end_tx,
            "fetch_start": min(s for s, e in exons),
            "fetch_end": max(e for s, e in exons)
        })

    return models


# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()

    offsets = load_offsets(args.offset)
    max_offset = max(offsets.values())
    bam = pysam.AlignmentFile(args.bam, "rb")

    models = load_tx_models(args.gtf, args.min_cds_len)
    start_tx = load_start_codons(args.gtf)
    models = [m for m in models if (m["gene"], m["tx"]) in start_tx]

    if not models:
        raise RuntimeError("No valid transcripts found")

    # --------------------------------------------------
    # PASS 1: select top transcripts (CDS only)
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

            rl = read.query_length
            if rl not in offsets:
                continue

            psite = (
                read.reference_start + offsets[rl]
                if strand == "+"
                else read.reference_end - offsets[rl] - 1
            )

            tx_pos = genomic_to_tx_spliced(psite, m["exons"], strand)
            if tx_pos is None:
                continue
            if m["cds_start_tx"] <= tx_pos < m["cds_end_tx"]:
                tx_counts[tx_id] += 1

    top_tx = {
        tx for tx, _ in sorted(
            tx_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:args.top_tx]
    }

    # --------------------------------------------------
    # PASS 2: region-aware periodicity
    # --------------------------------------------------
    counts = defaultdict(int)

    for m in models:
        tx_id = (m["gene"], m["tx"])
        if tx_id not in top_tx:
            continue

        chrom, strand = m["chrom"], m["strand"]
        fs = max(0, m["fetch_start"] - max_offset)
        fe = m["fetch_end"] + max_offset

        for read in bam.fetch(chrom, fs, fe):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq:
                continue

            rl = read.query_length
            if rl not in offsets:
                continue

            psite = (
                read.reference_start + offsets[rl]
                if strand == "+"
                else read.reference_end - offsets[rl] - 1
            )

            tx_pos = genomic_to_tx_spliced(psite, m["exons"], strand)
            if tx_pos is None:
                continue

            rel = tx_pos - m["cds_start_tx"]
            frame = rel % 3

            if tx_pos < m["cds_start_tx"]:
                region = "5UTR"
            elif tx_pos < m["cds_end_tx"]:
                region = "CDS"
            else:
                region = "3UTR"

            counts[(region, rl, frame)] += 1

    # --------------------------------------------------
    # output table
    # --------------------------------------------------
    df = pd.DataFrame(
        [
            {
                "sample": args.sample,
                "region": region,
                "read_length": rl,
                "frame": frame,
                "count": c
            }
            for (region, rl, frame), c in counts.items()
        ]
    )

    out_tsv = f"{args.sample}.periodicity.by_region_by_length.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[periodicity] written {out_tsv}")

    if df.empty:
        return

    df["fraction"] = (
        df["count"] /
        df.groupby(["region", "read_length"])["count"].transform("sum")
    )

    # --------------------------------------------------
    # plotting: 3 heatmaps (5′UTR, CDS, 3′UTR)
    # --------------------------------------------------
    regions = ["5UTR", "CDS", "3UTR"]
    region_titles = ["5′ UTR", "CDS", "3′ UTR"]

    cmap = LinearSegmentedColormap.from_list("white_red", ["#ffffff", "#ff0000"])

    all_lengths = sorted(df["read_length"].unique())
    nrows = len(all_lengths)

    fig, axes = plt.subplots(
        ncols=3,
        figsize=(11, max(4, nrows * 0.35)),
        constrained_layout=True
    )

    im = None

    for i, (ax, region, title) in enumerate(zip(axes, regions, region_titles)):
        sub = df[df["region"] == region]

        pivot = (
            sub.pivot_table(
                index="read_length",
                columns="frame",
                values="fraction",
                fill_value=0.0
            )
            .reindex(all_lengths, fill_value=0.0)
            .sort_index()
        )

        im = ax.imshow(
            pivot.values,
            aspect="auto",
            cmap=cmap,
            vmin=0.0,
            vmax=1.0
        )

        # ---- titles ----
        ax.set_title(title)

        # ---- x axis: keep ticks, remove label ----
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["Frame 0", "Frame 1", "Frame 2"])
        ax.set_xlabel("")   # remove axis label

        # ---- y axis ----
        ax.set_yticks(range(nrows))
        ax.set_ylim(nrows - 0.5, -0.5)

        if i == 0:
            # left panel: show read length labels
            ax.yaxis.set_major_locator(FixedLocator(range(nrows)))
            ax.yaxis.set_major_formatter(
                FixedFormatter([str(x) for x in all_lengths])
            )
            ax.set_ylabel("Read length")
            ax.tick_params(axis="y", labelsize=11, pad=8)
        else:
            # CDS & 3′UTR: remove y ticks entirely
            ax.set_yticklabels([])
            ax.tick_params(axis="y", left=False)

    # ---- shared colorbar ----
    cbar = fig.colorbar(im, ax=axes, fraction=0.035, pad=0.02)
    cbar.set_label("Fraction")
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

    fig.suptitle(f"{args.sample} frame periodicity (CDS-anchored)")

    out_png = f"{args.sample}.periodicity.by_region_by_length.heatmap.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"[periodicity] written {out_png}")





if __name__ == "__main__":
    main()
