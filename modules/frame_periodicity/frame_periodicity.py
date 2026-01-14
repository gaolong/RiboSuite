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


# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Frame periodicity QC using CDS-aligned P-sites "
                    "(restricted to transcripts with annotated start codons)"
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
# GTF parsing helpers
# --------------------------------------------------
def parse_attrs(attr):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr))


def load_start_codons(gtf):
    """
    Return a set of (gene_id, transcript_id) with annotated start codons.
    """
    start_tx = set()

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            chrom, _, feature, start, end, _, strand, phase, attr = line.rstrip().split("\t")
            if feature != "start_codon":
                continue

            attrs = parse_attrs(attr)
            gene = attrs.get("gene_id")
            tx   = attrs.get("transcript_id")

            if gene and tx:
                start_tx.add((gene, tx))

    return start_tx


def load_tx_models(gtf, min_cds_len):
    """
    Load CDS models per transcript.
    Coordinates are stored in transcript order with cumulative CDS offsets.
    """
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

            if feature != "CDS":
                continue

            gene = attrs.get("gene_id")
            tx   = attrs.get("transcript_id")
            if not gene or not tx:
                continue

            if gene_type.get(gene) != "protein_coding":
                continue

            start0 = int(start) - 1
            end0   = int(end)
            ph = 0 if phase == "." else int(phase)

            tx_cds[(gene, tx, chrom, strand)].append((start0, end0, ph))

    models = []

    for (gene, tx, chrom, strand), exons in tx_cds.items():
        exons.sort(key=lambda x: x[0])

        cds_len = sum(e - s for s, e, _ in exons)
        if cds_len < min_cds_len:
            continue

        tx_pos = []
        cursor = 0
        for s, e, ph in exons:
            tx_pos.append((s, e, cursor, ph))
            cursor += (e - s)

        models.append({
            "gene": gene,
            "tx": tx,
            "chrom": chrom,
            "strand": strand,
            "exons": tx_pos,
            "cds_len": cds_len
        })

    return models


# --------------------------------------------------
# genomic → transcript coordinate
# --------------------------------------------------
def genomic_to_tx(psite, exons, strand):
    """
    Convert genomic P-site to CDS-relative transcript coordinate.
    Coordinate 0 == first nucleotide of annotated CDS start.
    """
    if strand == "+":
        for s, e, tx_start, _ in exons:
            if s <= psite < e:
                return tx_start + (psite - s)
    else:
        for s, e, tx_start, _ in reversed(exons):
            if s <= psite < e:
                return tx_start + (e - psite - 1)
    return None


# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()

    offsets = load_offsets(args.offset)
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Load transcript models
    models = load_tx_models(args.gtf, args.min_cds_len)

    # Restrict to transcripts with annotated start codons
    start_tx = load_start_codons(args.gtf)
    models = [
        m for m in models
        if (m["gene"], m["tx"]) in start_tx
    ]

    if not models:
        raise RuntimeError("No transcripts with annotated start codons found.")

    # --------------------------------------------------
    # PASS 1: select top transcripts by CDS P-site count
    # --------------------------------------------------
    tx_counts = defaultdict(int)

    for m in models:
        chrom, strand = m["chrom"], m["strand"]
        tx_id = (m["gene"], m["tx"])

        for read in bam.fetch(chrom):
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

            tx_pos = genomic_to_tx(psite, m["exons"], strand)
            if tx_pos is not None:
                tx_counts[tx_id] += 1

    top_tx = {
        tx for tx, _ in sorted(
            tx_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:args.top_tx]
    }

    # --------------------------------------------------
    # PASS 2: frame periodicity (anchored at CDS start)
    # --------------------------------------------------
    counts = defaultdict(int)

    for m in models:
        tx_id = (m["gene"], m["tx"])
        if tx_id not in top_tx:
            continue

        chrom, strand = m["chrom"], m["strand"]

        for read in bam.fetch(chrom):
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

            tx_pos = genomic_to_tx(psite, m["exons"], strand)
            if tx_pos is None:
                continue

            frame = tx_pos % 3
            counts[(rl, frame)] += 1

    # --------------------------------------------------
    # output table
    # --------------------------------------------------
    rows = [
        {
            "sample": args.sample,
            "read_length": rl,
            "frame": fr,
            "count": c
        }
        for (rl, fr), c in counts.items()
    ]

    df = pd.DataFrame(rows)

    out_tsv = f"{args.sample}.periodicity.by_region_by_length.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[periodicity] written {out_tsv}")

    if df.empty:
        print("[periodicity] no data available for plotting")
        return

    df["fraction"] = (
        df["count"] /
        df.groupby("read_length")["count"].transform("sum")
    )

    # --------------------------------------------------
    # heatmap
    # --------------------------------------------------
    pivot = (
        df.pivot_table(
            index="read_length",
            columns="frame",
            values="fraction",
            fill_value=0.0
        )
        .sort_index()
    )

    fig, ax = plt.subplots(figsize=(4, max(3, len(pivot) * 0.25)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis")

    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Frame 0", "Frame 1", "Frame 2"])

    ax.set_xlabel("Frame")
    ax.set_ylabel("Read length")
    ax.set_title(f"{args.sample} frame periodicity (annotated starts)")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Fraction")

    out_png = f"{args.sample}.periodicity.by_region_by_length.heatmap.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"[periodicity] written {out_png}")


if __name__ == "__main__":
    main()
