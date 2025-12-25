import os
import pysam
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--sample", required=True)
    return ap.parse_args()

def load_offsets(path):
    df = pd.read_csv(path, sep="\t")

    if "psite_offset" in df.columns:
        offset_col = "psite_offset"
    elif "offset" in df.columns:
        offset_col = "offset"
    else:
        raise ValueError(
            f"Offset column not found in {path}. "
            f"Columns: {list(df.columns)}"
        )

    return dict(zip(df["read_length"], df[offset_col]))

def load_cds(gtf):
    cds = defaultdict(list)
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "CDS":
                continue
            chrom, _, _, start, end, _, strand, frame, attr = fields
            start, end = int(start)-1, int(end)
            cds[(chrom, strand)].append((start, end))
    return cds

def assign_frame(psite, cds_intervals):
    for start, end in cds_intervals:
        if start <= psite < end:
            return (psite - start) % 3
    return None

def main():
    args = parse_args()

    offsets = load_offsets(args.offset)
    cds_map = load_cds(args.gtf)

    bam = pysam.AlignmentFile(args.bam, "rb")

    frame_counts = defaultdict(int)

    for read in bam.fetch():
        if read.is_unmapped:
            continue

        read_len = read.query_length
        if read_len not in offsets:
            continue

        offset = offsets[read_len]
        strand = "-" if read.is_reverse else "+"
        chrom = bam.get_reference_name(read.reference_id)

        psite = read.reference_start + offset if strand == "+" else read.reference_end - offset - 1

        key = (chrom, strand)
        if key not in cds_map:
            continue

        frame = assign_frame(psite, cds_map[key])
        if frame is not None:
            frame_counts[frame] += 1

    # ---- output tables ----
    df = pd.DataFrame(
        [{"frame": k, "count": v} for k, v in frame_counts.items()]
    ).sort_values("frame")

    df.to_csv(f"{args.sample}.frame_counts.tsv", sep="\t", index=False)

    df["fraction"] = df["count"] / df["count"].sum()
    df[["frame", "fraction"]].to_csv(
        f"{args.sample}.frame_fraction.tsv", sep="\t", index=False
    )

    # ---- plot ----
    plt.bar(df["frame"], df["fraction"])
    plt.xlabel("Reading frame")
    plt.ylabel("Fraction of P-sites")
    plt.title(f"{args.sample} frame periodicity")
    plt.savefig(f"{args.sample}.frame_periodicity.png", dpi=150, bbox_inches="tight")

if __name__ == "__main__":
    main()
