#!/usr/bin/env python3

import pysam
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import re

# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Frame periodicity QC using top-N protein-coding transcripts by CDS P-site coverage"
    )
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min_cds_len", type=int, default=150)
    ap.add_argument("--top_tx", type=int, default=1000,
                    help="Keep top N transcripts by CDS P-site count")
    return ap.parse_args()

# --------------------------------------------------
# offsets
# --------------------------------------------------
def load_offsets(path):
    df = pd.read_csv(path, sep="\t")
    col = "psite_offset" if "psite_offset" in df.columns else "offset"
    return dict(zip(df["read_length"], df[col]))

# --------------------------------------------------
# GTF parsing
# --------------------------------------------------
def parse_attrs(attr):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr))

def load_all_tx_models(gtf, min_cds_len):
    tx_cds = defaultdict(list)
    gene_type = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, phase, attr = line.rstrip().split("\t")
            attrs = parse_attrs(attr)

            if feature == "gene":
                gene_type[attrs["gene_id"]] = attrs.get("gene_type") or attrs.get("gene_biotype")

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

        cds_start = exons[0][0] if strand == "+" else exons[-1][1] - 1
        cds_end   = exons[-1][1] - 1 if strand == "+" else exons[0][0]

        models.append({
            "gene": gene,
            "tx": tx,
            "chrom": chrom,
            "strand": strand,
            "cds_exons": exons,
            "cds_start": cds_start,
            "cds_end": cds_end
        })

    return models

# --------------------------------------------------
# utilities
# --------------------------------------------------
def frame_from_cds_start(psite, cds_start, strand):
    return (psite - cds_start) % 3 if strand == "+" else (cds_start - psite) % 3

# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()
    offsets = load_offsets(args.offset)
    bam = pysam.AlignmentFile(args.bam, "rb")

    models = load_all_tx_models(args.gtf, args.min_cds_len)

    # --------------------------------------------------
    # PASS 1: count CDS P-sites per transcript
    # --------------------------------------------------
    tx_psite_counts = defaultdict(int)

    for m in models:
        chrom, strand = m["chrom"], m["strand"]
        tx_id = (m["gene"], m["tx"])

        for s, e, _ in m["cds_exons"]:
            for read in bam.fetch(chrom, s, e):
                if read.is_unmapped:
                    continue
                rl = read.query_length
                if rl not in offsets:
                    continue

                psite = (
                    read.reference_start + offsets[rl]
                    if strand == "+"
                    else read.reference_end - offsets[rl] - 1
                )
                if s <= psite < e:
                    tx_psite_counts[tx_id] += 1

    # select top N transcripts
    top_tx = set(
        tx for tx, _ in sorted(
            tx_psite_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:args.top_tx]
    )

    # --------------------------------------------------
    # PASS 2: frame periodicity on top transcripts
    # --------------------------------------------------
    counts = defaultdict(int)

    for m in models:
        tx_id = (m["gene"], m["tx"])
        if tx_id not in top_tx:
            continue

        chrom, strand = m["chrom"], m["strand"]

        for s, e, _ in m["cds_exons"]:
            for read in bam.fetch(chrom, s, e):
                if read.is_unmapped:
                    continue
                rl = read.query_length
                if rl not in offsets:
                    continue

                psite = (
                    read.reference_start + offsets[rl]
                    if strand == "+"
                    else read.reference_end - offsets[rl] - 1
                )
                if not (s <= psite < e):
                    continue

                frame = frame_from_cds_start(psite, m["cds_start"], strand)
                counts[("CDS", rl, frame)] += 1

    # --------------------------------------------------
    # table
    # --------------------------------------------------
    rows = [
        dict(
            sample=args.sample,
            region="CDS",
            read_length=rl,
            frame=fr,
            count=c
        )
        for (region, rl, fr), c in counts.items()
    ]

    df = pd.DataFrame(rows)
    df["fraction"] = df["count"] / df.groupby(
        ["region", "read_length"]
    )["count"].transform("sum")

    out_tsv = f"{args.sample}.periodicity.by_region_by_length.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    # --------------------------------------------------
    # heatmap
    # --------------------------------------------------
    fig, ax = plt.subplots(figsize=(4, 5))

    mat = (
        df.pivot(index="read_length", columns="frame", values="fraction")
        .fillna(0)
        .sort_index()
    )

    vmax = float(mat.to_numpy().max()) if mat.size else 0.0
    if vmax <= 0:
        vmax = 1.0  # fallback to avoid a zero-range colorbar
    im = ax.imshow(mat, aspect="auto", origin="lower",
                   cmap="Reds", vmin=0, vmax=vmax)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Frame 0", "Frame 1", "Frame 2"])
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index)
    ax.set_title("CDS (top transcripts)")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("P-site fraction", rotation=270, labelpad=12)

    fig.savefig(
        f"{args.sample}.periodicity.by_region_by_length.heatmap.png",
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()


if __name__ == "__main__":
    main()
