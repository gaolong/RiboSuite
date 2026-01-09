#!/usr/bin/env python3

import pysam
import argparse
import pandas as pd
import sys
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert BAM to 1-nt P-site BED files split by read length"
    )
    p.add_argument("--bam", required=True)
    p.add_argument("--offset", required=True,
                   help="TSV with columns: read_length, psite_offset")
    p.add_argument("--out", required=True,
                   help="Output prefix (e.g. sample.psite)")
    p.add_argument("--mapq", type=int, default=100)
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
out_handles = {}
written_by_len = defaultdict(int)

n_unmapped = 0
n_low_mapq = 0
n_no_offset = 0
n_written_total = 0

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
    else:
        psite = r.reference_end - 1 - offset

    if psite < 0:
        continue

    # Lazily open output file for this read length
    if L not in out_handles:
        out_path = f"{args.out}.len{L}.bed"
        out_handles[L] = open(out_path, "w")

    out_handles[L].write(
        f"{bam.get_reference_name(r.reference_id)}\t"
        f"{psite}\t{psite + 1}\n"
    )

    written_by_len[L] += 1
    n_written_total += 1

# --------------------------------------------------
# Close files
# --------------------------------------------------
for fh in out_handles.values():
    fh.close()

# --------------------------------------------------
# Summary (stderr)
# --------------------------------------------------
print(
    f"[bam_to_psite_by_len] "
    f"unmapped/secondary={n_unmapped}, "
    f"low_mapq={n_low_mapq}, "
    f"no_offset={n_no_offset}, "
    f"written_total={n_written_total}",
    file=sys.stderr
)

for L in sorted(written_by_len):
    print(
        f"[bam_to_psite_by_len] "
        f"len={L}\twritten={written_by_len[L]}",
        file=sys.stderr
    )
