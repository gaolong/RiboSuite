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
        description="Frame periodicity by region (5UTR/CDS/3UTR) with CDS start-window logic preserved"
    )
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--sample", required=True)

    ap.add_argument("--start_window", default="-30,10",
                    help="Start window relative to CDS start (nt), e.g. -30,10")
    ap.add_argument("--min_cds_len", type=int, default=150,
                    help="Minimum CDS length to keep transcript")

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

def load_longest_tx_models(gtf, min_cds_len):
    """
    One protein-coding transcript (longest CDS) per gene.
    """
    tx_cds = defaultdict(list)
    tx_bounds = {}
    gene_type = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, phase, attr = line.rstrip().split("\t")
            attrs = parse_attrs(attr)

            if feature == "gene":
                gene_type[attrs["gene_id"]] = attrs.get("gene_type") or attrs.get("gene_biotype")

            if feature not in ("exon", "CDS"):
                continue

            gene = attrs.get("gene_id")
            tx = attrs.get("transcript_id")
            if not gene or not tx:
                continue

            start0 = int(start) - 1
            end0 = int(end)

            tx_bounds.setdefault((gene, tx, chrom, strand), [start0, end0])
            tx_bounds[(gene, tx, chrom, strand)][0] = min(tx_bounds[(gene, tx, chrom, strand)][0], start0)
            tx_bounds[(gene, tx, chrom, strand)][1] = max(tx_bounds[(gene, tx, chrom, strand)][1], end0)

            if feature == "CDS":
                ph = 0 if phase == "." else int(phase)
                tx_cds[(gene, tx, chrom, strand)].append((start0, end0, ph))

    best = {}

    for (gene, tx, chrom, strand), exons in tx_cds.items():
        if gene_type.get(gene) != "protein_coding":
            continue

        exons.sort(key=lambda x: x[0])
        cds_len = sum(e - s for s, e, _ in exons)
        if cds_len < min_cds_len:
            continue

        if gene not in best or cds_len > best[gene]["cds_len"]:
            cds_start = exons[0][0] if strand == "+" else exons[-1][1] - 1
            cds_end = exons[-1][1] - 1 if strand == "+" else exons[0][0]
            tx_start, tx_end = tx_bounds[(gene, tx, chrom, strand)]

            best[gene] = {
                "chrom": chrom,
                "strand": strand,
                "cds_exons": exons,
                "cds_start": cds_start,
                "cds_end": cds_end,
                "tx_start": tx_start,
                "tx_end": tx_end,
                "cds_len": cds_len
            }

    return list(best.values())

# --------------------------------------------------
# utilities
# --------------------------------------------------
def frame_from_cds_start(psite, cds_start, strand):
    return (psite - cds_start) % 3 if strand == "+" else (cds_start - psite) % 3

def build_start_window_segments(cds_exons, cds_start, strand, wL, wR):
    segs = []
    if strand == "+":
        cursor = cds_start + wL
        remaining = wR - wL + 1
        for s, e, ph in cds_exons:
            if cursor >= e:
                continue
            if cursor < s:
                cursor = s
            take = min(e - cursor, remaining)
            if take > 0:
                segs.append((cursor, cursor + take, ph))
                remaining -= take
                cursor += take
            if remaining <= 0:
                break
    else:
        cursor = cds_start - wL
        remaining = wR - wL + 1
        for s, e, ph in reversed(cds_exons):
            if cursor < s:
                continue
            if cursor >= e:
                cursor = e - 1
            take = min(cursor - s + 1, remaining)
            if take > 0:
                segs.append((cursor - take + 1, cursor + 1, ph))
                remaining -= take
                cursor -= take
            if remaining <= 0:
                break
    return segs

# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()
    offsets = load_offsets(args.offset)
    wL, wR = map(int, args.start_window.split(","))

    models = load_longest_tx_models(args.gtf, args.min_cds_len)
    bam = pysam.AlignmentFile(args.bam, "rb")

    counts = defaultdict(int)

    for m in models:
        chrom, strand = m["chrom"], m["strand"]

        # ---------- CDS (START WINDOW ONLY, UNCHANGED LOGIC) ----------
        cds_segs = build_start_window_segments(
            m["cds_exons"], m["cds_start"], strand, wL, wR
        )

        for s, e, ph in cds_segs:
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

        # ---------- 5UTR ----------
        utr5 = (
            (m["tx_start"], m["cds_start"])
            if strand == "+"
            else (m["cds_start"] + 1, m["tx_end"])
        )
        if utr5[0] < utr5[1]:
            for read in bam.fetch(chrom, utr5[0], utr5[1]):
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
                if not (utr5[0] <= psite < utr5[1]):
                    continue
                frame = frame_from_cds_start(psite, m["cds_start"], strand)
                counts[("5UTR", rl, frame)] += 1

        # ---------- 3UTR ----------
        utr3 = (
            (m["cds_end"] + 1, m["tx_end"])
            if strand == "+"
            else (m["tx_start"], m["cds_end"])
        )
        if utr3[0] < utr3[1]:
            for read in bam.fetch(chrom, utr3[0], utr3[1]):
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
                if not (utr3[0] <= psite < utr3[1]):
                    continue
                frame = frame_from_cds_start(psite, m["cds_start"], strand)
                counts[("3UTR", rl, frame)] += 1

    # --------------------------------------------------
    # table
    # --------------------------------------------------
    rows = [
        {
            "sample": args.sample,
            "region": region,
            "read_length": rl,
            "frame": fr,
            "count": c
        }
        for (region, rl, fr), c in counts.items()
    ]

    df = pd.DataFrame(rows)
    df["fraction"] = df["count"] / df.groupby(
        ["region", "read_length"]
    )["count"].transform("sum")

    out_tsv = f"{args.sample}.periodicity.by_region_by_length.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    # --------------------------------------------------
    # heatmaps (3 panels, frame 0/1/2)
    # --------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    regions = ["5UTR", "CDS", "3UTR"]

    title_map = {
        "5UTR": "5' UTR",
        "CDS": "CDS",
        "3UTR": "3' UTR"
    }

    for ax, region in zip(axes, regions):
        sub = df[df["region"] == region]
        mat = (
            sub.pivot(index="read_length", columns="frame", values="fraction")
            .fillna(0)
            .sort_index()
        )
        im = ax.imshow(
            mat,
            aspect="auto",
            origin="lower",
            cmap="Reds",
            vmin=0,
            vmax=1
        )
        ax.set_title(title_map[region])
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["Frame 0", "Frame 1", "Frame 2"])
        ax.set_yticks(range(len(mat.index)))
        ax.set_yticklabels(mat.index)

    # shared y-axis label
    fig.text(0.04, 0.5, "Read length", va="center", rotation="vertical")

    # dedicated colorbar axis
    # dedicated colorbar axis
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("P-site fraction", rotation=270, labelpad=12)
    cbar.ax.yaxis.label.set_verticalalignment("center")

    fig.suptitle(f"{args.sample} frame periodicity by region")
    plt.subplots_adjust(wspace=0.25, right=0.9)
    fig.savefig(
        f"{args.sample}.periodicity.by_region_by_length.heatmap.png",
        dpi=150
    )
    plt.close()



if __name__ == "__main__":
    main()
