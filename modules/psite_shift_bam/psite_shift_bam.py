#!/usr/bin/env python3

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
import pysam

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../../"))
BIN_DIR = os.path.join(ROOT_DIR, "bin")
sys.path.insert(0, BIN_DIR)

from psite_utils import load_offsets, get_psite


def parse_args():
    p = argparse.ArgumentParser(
        description="Rewrite BAM into P-site-shifted BAM so downstream zero-offset logic is valid."
    )
    p.add_argument("--bam", required=True, help="Input BAM")
    p.add_argument("--offsets", required=True, help="Input per-length P-site offsets TSV")
    p.add_argument("--sample", required=True, help="Sample ID")
    p.add_argument("--out_bam", required=True, help="Output shifted BAM")
    p.add_argument("--out_offsets", required=True, help="Output zero-offset TSV")
    return p.parse_args()


def make_zero_offset_table(read_lengths, sample_id, out_path):
    df = pd.DataFrame({
        "sample_id": sample_id,
        "read_length": sorted(read_lengths),
        "psite_offset": 0,
    })
    df.to_csv(out_path, sep="\t", index=False)


def main():
    args = parse_args()

    offsets = load_offsets(args.offsets)
    bam_in = pysam.AlignmentFile(args.bam, "rb")

    header = bam_in.header.to_dict()
    bam_out = pysam.AlignmentFile(args.out_bam, "wb", header=header)

    used_lengths = set()
    stats = defaultdict(int)

    for read in bam_in.fetch(until_eof=True):
        stats["reads_seen"] += 1

        if read.is_unmapped:
            stats["unmapped"] += 1
            continue
        if read.is_secondary or read.is_supplementary:
            stats["secondary_or_supp"] += 1
            continue

        psite, read_len = get_psite(read, offsets)
        if psite is None or read_len is None:
            stats["no_psite"] += 1
            continue

        used_lengths.add(int(read_len))

        new_read = pysam.AlignedSegment.from_dict(read.to_dict(), header=bam_out.header)

        # Keep original read length and orientation.
        # Goal:
        #   plus strand: reference_start == psite
        #   minus strand: reference_end - 1 == psite
        #
        # We encode the alignment as a simple contiguous match of read_len M.
        # This is a simplified representation intended for downstream P-site-based analyses.
        new_read.cigartuples = [(0, int(read_len))]  # M

        if read.is_reverse:
            new_start = int(psite) - int(read_len) + 1
        else:
            new_start = int(psite)

        if new_start < 0:
            stats["negative_start_skipped"] += 1
            continue

        new_read.reference_start = new_start

        # Remove tags that may become inconsistent after rewriting
        for tag in ("MD", "NM"):
            try:
                new_read.set_tag(tag, None)
            except Exception:
                pass

        bam_out.write(new_read)
        stats["written"] += 1

    bam_in.close()
    bam_out.close()

    if not used_lengths:
        raise RuntimeError("No reads written to shifted BAM; cannot build zero-offset table")

    make_zero_offset_table(used_lengths, args.sample, args.out_offsets)

    sys.stderr.write(
        "[psite_shift_bam] "
        + " ".join(f"{k}={v}" for k, v in sorted(stats.items()))
        + "\n"
    )


if __name__ == "__main__":
    main()