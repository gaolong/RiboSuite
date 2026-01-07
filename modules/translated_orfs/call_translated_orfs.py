#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
from collections import defaultdict, Counter
from pathlib import Path

# -------------------------
# CLI
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Rule-based identification of translated ORFs from Ribo-seq P-site signal"
    )
    p.add_argument("--bam", required=True, help="Deduplicated, coordinate-sorted BAM")
    p.add_argument("--gtf", required=True, help="Annotation GTF")
    p.add_argument("--psite_offsets", required=True, help="P-site offset table")
    p.add_argument("--out_prefix", required=True)

    # Conservative defaults
    p.add_argument("--min_psites", type=int, default=30)
    p.add_argument("--min_len_codons", type=int, default=30)
    p.add_argument("--min_frame_fraction", type=float, default=0.7)
    p.add_argument("--allow_non_aug", action="store_true")

    return p.parse_args()

# -------------------------
# Loaders
# -------------------------
def load_psite_offsets(path):
    """
    Expect columns:
      read_length, psite_offset
    """
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["read_length"], df["psite_offset"]))


def load_cds_from_gtf(gtf_path):
    """
    v1: extract annotated CDS per transcript
    Returns:
      dict transcript_id -> list of (chrom, start, end, strand, frame)
    """
    cds = defaultdict(list)

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "CDS":
                continue

            chrom, _, _, start, end, _, strand, frame, attrs = fields
            start, end = int(start), int(end)

            attr = dict(
                item.strip().replace('"', '').split(" ")
                for item in attrs.strip(";").split("; ")
            )
            tid = attr.get("transcript_id")
            gid = attr.get("gene_id")

            cds[tid].append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "frame": frame,
                "gene_id": gid
            })

    return cds

# -------------------------
# Core logic
# -------------------------
def iter_psites(bam, offsets):
    """
    Yield (chrom, psite_pos, strand, read_length)
    """
    for read in bam.fetch():
        if read.is_unmapped:
            continue
        if read.query_length not in offsets:
            continue

        offset = offsets[read.query_length]
        strand = "-" if read.is_reverse else "+"

        if strand == "+":
            psite = read.reference_start + offset
        else:
            psite = read.reference_end - offset - 1

        yield read.reference_name, psite, strand, read.query_length


def assign_psites_to_orf(psites, cds_blocks):
    """
    Count P-sites within CDS blocks and compute frame statistics.
    """
    frame_counts = Counter()
    total = 0

    for chrom, pos, strand in psites:
        for block in cds_blocks:
            if block["chrom"] != chrom or block["strand"] != strand:
                continue
            if block["start"] <= pos <= block["end"]:
                frame = (pos - block["start"]) % 3
                frame_counts[frame] += 1
                total += 1

    return total, frame_counts

# -------------------------
# Rules
# -------------------------
def rule_min_psites(total, min_psites):
    return total >= min_psites


def rule_frame_dominance(frame_counts, min_fraction):
    if not frame_counts:
        return False, 0.0
    top = max(frame_counts.values())
    frac = top / sum(frame_counts.values())
    return frac >= min_fraction, frac


# -------------------------
# Main
# -------------------------
def main():
    args = parse_args()

    offsets = load_psite_offsets(args.psite_offsets)
    cds_by_tx = load_cds_from_gtf(args.gtf)
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Precompute P-sites
    psites = list(iter_psites(bam, offsets))

    rows = []

    for tx, blocks in cds_by_tx.items():
        gene_id = blocks[0]["gene_id"]
        total, frame_counts = assign_psites_to_orf(psites, blocks)

        rules_passed = []
        rules_failed = []

        if rule_min_psites(total, args.min_psites):
            rules_passed.append("min_psites")
        else:
            rules_failed.append("min_psites")

        ok_frame, frac = rule_frame_dominance(
            frame_counts, args.min_frame_fraction
        )
        if ok_frame:
            rules_passed.append("frame_dominance")
        else:
            rules_failed.append("frame_dominance")

        rows.append({
            "gene_id": gene_id,
            "transcript_id": tx,
            "psite_count": total,
            "frame_fraction": round(frac, 3),
            "rules_passed": ",".join(rules_passed),
            "rules_failed": ",".join(rules_failed)
        })

    df = pd.DataFrame(rows)
    df.to_csv(f"{args.out_prefix}.translated_orfs.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
