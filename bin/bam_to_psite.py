#!/usr/bin/env python3

import pysam
import argparse
import pandas as pd
import sys
from collections import defaultdict

# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Convert BAM to 1-nt P-site BED files. "
            "Always produces a combined BED; optionally split by read length; "
            "optionally generate strand-specific BEDs."
        )
    )
    p.add_argument("--bam", required=True)
    p.add_argument(
        "--offset",
        required=True,
        help="TSV with columns: read_length, psite_offset"
    )
    p.add_argument(
        "--out",
        required=True,
        help="Output prefix (e.g. sample.psite)"
    )
    p.add_argument("--mapq", type=int, default=100)
    p.add_argument(
        "--by_length",
        action="store_true",
        help="Also generate per-read-length BED files (*.len{L}.bed)"
    )
    p.add_argument(
        "--strand_specific",
        action="store_true",
        help="Also generate strand-specific BED files (*.pos.bed, *.neg.bed). "
             "pos = reads on '+' strand (not reverse); neg = reads on '-' strand (reverse)."
    )
    return p.parse_args()

args = parse_args()

# --------------------------------------------------
# Load offsets: read_length -> psite_offset
# --------------------------------------------------
off = pd.read_csv(args.offset, sep="\t")

required_cols = {"read_length", "psite_offset"}
if not required_cols.issubset(off.columns):
    sys.exit(
        f"ERROR: offset file must contain columns {required_cols}, "
        f"found {set(off.columns)}"
    )

off["read_length"] = off["read_length"].astype(int)
off["psite_offset"] = off["psite_offset"].astype(int)
offsets = dict(zip(off["read_length"], off["psite_offset"]))

# --------------------------------------------------
# Open BAM
# --------------------------------------------------
bam = pysam.AlignmentFile(args.bam, "rb")

# --------------------------------------------------
# Output handles & counters
# --------------------------------------------------
out_handles = {}                  # per-length BEDs (combined)
out_handles_pos = {}              # per-length BEDs (pos strand) if enabled
out_handles_neg = {}              # per-length BEDs (neg strand) if enabled

written_by_len = defaultdict(int)
written_by_len_pos = defaultdict(int)
written_by_len_neg = defaultdict(int)

# Combined BED (ALWAYS written)
out_all = open(f"{args.out}.all.bed", "w")

# Strand-specific BEDs (optional)
out_pos = open(f"{args.out}.pos.bed", "w") if args.strand_specific else None
out_neg = open(f"{args.out}.neg.bed", "w") if args.strand_specific else None

n_unmapped = 0
n_low_mapq = 0
n_no_offset = 0
n_written_total = 0
n_written_pos = 0
n_written_neg = 0

# --------------------------------------------------
# Main loop
# --------------------------------------------------
for r in bam.fetch(until_eof=True):

    # Skip unmapped / secondary / supplementary
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        n_unmapped += 1
        continue

    # MAPQ filter
    if r.mapping_quality < args.mapq:
        n_low_mapq += 1
        continue

    # Length-specific offset
    L = r.query_length
    if L not in offsets:
        n_no_offset += 1
        continue

    offset = offsets[L]

    # Compute P-site from 5' end
    if not r.is_reverse:
        psite = r.reference_start + offset
        is_pos = True
    else:
        psite = r.reference_end - 1 - offset
        is_pos = False

    if psite < 0:
        continue

    chrom = bam.get_reference_name(r.reference_id)
    bed_line = f"{chrom}\t{psite}\t{psite + 1}\n"

    # Combined BED (always)
    out_all.write(bed_line)

    # Optional strand-specific BEDs
    if args.strand_specific:
        if is_pos:
            out_pos.write(bed_line)
            n_written_pos += 1
        else:
            out_neg.write(bed_line)
            n_written_neg += 1

    # Optional per-length BEDs (combined and optionally strand-specific)
    if args.by_length:
        # combined
        if L not in out_handles:
            out_handles[L] = open(f"{args.out}.len{L}.bed", "w")
        out_handles[L].write(bed_line)
        written_by_len[L] += 1

        # strand-specific per length
        if args.strand_specific:
            if is_pos:
                if L not in out_handles_pos:
                    out_handles_pos[L] = open(f"{args.out}.pos.len{L}.bed", "w")
                out_handles_pos[L].write(bed_line)
                written_by_len_pos[L] += 1
            else:
                if L not in out_handles_neg:
                    out_handles_neg[L] = open(f"{args.out}.neg.len{L}.bed", "w")
                out_handles_neg[L].write(bed_line)
                written_by_len_neg[L] += 1

    n_written_total += 1

# --------------------------------------------------
# Close files
# --------------------------------------------------
for fh in out_handles.values():
    fh.close()
for fh in out_handles_pos.values():
    fh.close()
for fh in out_handles_neg.values():
    fh.close()

out_all.close()
if out_pos:
    out_pos.close()
if out_neg:
    out_neg.close()

# --------------------------------------------------
# Summary (stderr)
# --------------------------------------------------
print(
    f"[bam_to_psite] "
    f"unmapped/secondary={n_unmapped}, "
    f"low_mapq={n_low_mapq}, "
    f"no_offset={n_no_offset}, "
    f"written_total={n_written_total}"
    + (f", written_pos={n_written_pos}, written_neg={n_written_neg}" if args.strand_specific else ""),
    file=sys.stderr
)

if args.by_length:
    for L in sorted(written_by_len):
        msg = f"[bam_to_psite] len={L}\twritten={written_by_len[L]}"
        if args.strand_specific:
            msg += f"\tpos={written_by_len_pos[L]}\tneg={written_by_len_neg[L]}"
        print(msg, file=sys.stderr)