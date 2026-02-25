#!/usr/bin/env python3
"""
CDS-level P-site quantification (featureCounts-aligned semantics)

Goal: produce gene-level CDS counts from Ribo-seq BAM that are compatible with
RNA-seq gene counts from featureCounts for TE analysis.

Counting contract (aligned to featureCounts defaults, simplified):
- Count each read once (primary alignment only; optional MAPQ filter)
- Assign by overlap of the read's P-site position with gene CDS union regions
- Strand-consistent assignment (read orientation must match gene strand)
- If P-site overlaps multiple genes: skip (default) or count all with --skip_multi_overlap
- Optional --inframe_only filtering relative to CDS start (not GTF exon frame)
"""

import argparse
import os
import re
import sys
from bisect import bisect_right
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import pandas as pd
import pysam

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../../"))
BIN_DIR = os.path.join(ROOT_DIR, "bin")
sys.path.insert(0, BIN_DIR)

from psite_utils import get_psite  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(
        description="Gene-level CDS P-site quantification (featureCounts-aligned)"
    )
    p.add_argument("--bam", required=True, help="STAR-aligned BAM")
    p.add_argument("--gtf", required=True, help="GTF annotation")

    # renamed flag
    p.add_argument(
        "--psite_offsets",
        required=True,
        help="P-site offsets TSV (must contain read_length and psite_offset; optional sample_id)",
    )

    p.add_argument("--sample", required=True, help="Sample ID")
    p.add_argument("--out_prefix", default="cds_quant", help="Output prefix")

    # NOTE: STAR MAPQ can be 255 for unique; default 100 is conservative.
    p.add_argument("--min_mapq", type=int, default=100, help="Minimum MAPQ (default: 100)")

    p.add_argument(
        "--include_duplicates",
        action="store_true",
        help="Include PCR/optical duplicates (default: excluded if flagged)",
    )
    p.add_argument(
        "--skip_multi_overlap",
        action="store_true",
        help="If P-site overlaps multiple genes, skip it as ambiguous (default: count all)",
    )
    p.add_argument(
        "--inframe_only",
        action="store_true",
        help="Only count frame-0 P-sites relative to CDS start (recommended for ribo QC/TE)",
    )
    return p.parse_args()


def load_offsets(offset_file: str) -> Dict[int, int]:
    """
    Load P-site offsets from a TSV with columns:
      read_length, psite_offset
    Optional:
      sample_id (if multiple, uses the first sample_id and warns)

    Returns:
      dict: {read_length (int) -> psite_offset (int)}
    """
    df = pd.read_csv(offset_file, sep="\t")

    required = {"read_length", "psite_offset"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Offsets file must contain columns {required}. "
            f"Found: {list(df.columns)}"
        )

    if "sample_id" in df.columns:
        sample_ids = df["sample_id"].dropna().unique()
        if len(sample_ids) > 1:
            print(
                f"[WARN] Multiple sample_id values in offsets file: {sample_ids}. "
                f"Using offsets from sample_id={sample_ids[0]}",
                file=sys.stderr,
            )
        if len(sample_ids) >= 1:
            df = df[df["sample_id"] == sample_ids[0]]

    offsets = dict(
        zip(
            df["read_length"].astype(int),
            df["psite_offset"].astype(int),
        )
    )
    if not offsets:
        raise ValueError("No valid offsets loaded from offsets file")
    return offsets


def parse_gtf_cds_union(gtf_file: str):
    """
    Parse CDS entries from a GTF and build gene-level CDS union intervals.

    Returns:
      gene_meta: dict[gene_id] = {
          "chrom": str,
          "strand": str,
          "cds_start": int,
          "cds_len": int,
          "gene_name": str,
          "intervals": [(s,e),...]
      }
      chrom_index: dict[chrom] = {
            "starts": [start0, start1, ...]  (sorted ints),
            "entries": [(start, end, gene_id, strand), ...] (same order),
            "max_len": int (max interval length on chrom)
      }

    Notes:
    - Coordinates stored as 0-based start, end-exclusive.
    - cds_start is the genomic coordinate of the first CDS base of the ORF:
        '+' : minimum CDS start
        '-' : maximum CDS end - 1
    """
    attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
    cds_raw: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = defaultdict(list)
    gene_names: Dict[str, str] = {}

    with open(gtf_file, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attr = parts
            if feature != "CDS":
                continue

            attrs = dict(attr_re.findall(attr))
            gene = attrs.get("gene_id")
            if not gene:
                continue

            gene_name = attrs.get("gene_name")
            if gene_name and gene not in gene_names:
                gene_names[gene] = gene_name

            # GTF is 1-based inclusive; convert to 0-based half-open
            s = int(start) - 1
            e = int(end)
            if e <= s:
                continue

            cds_raw[(gene, chrom, strand)].append((s, e))

    gene_meta = {}
    for (gene, chrom, strand), intervals in cds_raw.items():
        intervals.sort()
        merged: List[Tuple[int, int]] = []
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e:  # overlap/adjacent
                cur_e = max(cur_e, e)
            else:
                merged.append((cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((cur_s, cur_e))

        if strand == "+":
            cds_start = min(s for s, _ in merged)
        else:
            cds_start = max(e for _, e in merged) - 1

        cds_len = sum(e - s for s, e in merged)
        gene_meta[gene] = {
            "chrom": chrom,
            "strand": strand,
            "cds_start": int(cds_start),
            "cds_len": int(cds_len),
            "gene_name": gene_names.get(gene, ""),
            "intervals": merged,
        }

    # Build a per-chrom interval index for fast overlap queries by point (psite)
    # Defensive: ensure starts is always a list[int], never tuples.
    chrom_entries: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)
    for gene, meta in gene_meta.items():
        chrom = meta["chrom"]
        strand = meta["strand"]
        for iv in meta["intervals"]:
            if not (isinstance(iv, (tuple, list)) and len(iv) == 2):
                raise ValueError(f"Bad interval for gene {gene}: {iv} (expected (s,e))")
            s, e = iv
            s = int(s)
            e = int(e)
            chrom_entries[chrom].append((s, e, gene, strand))

    chrom_index = {}
    for chrom, entries in chrom_entries.items():
        bad = [
            x for x in entries
            if not (isinstance(x, tuple) and len(x) == 4 and isinstance(x[0], int) and isinstance(x[1], int))
        ]
        if bad:
            raise ValueError(f"Bad chrom_entries on {chrom}. Example: {bad[0]}")

        entries.sort(key=lambda x: x[0])
        starts = [s for (s, e, gene, strand) in entries]  # guaranteed ints

        max_len = 0
        for s, e, _, _ in entries:
            max_len = max(max_len, e - s)

        chrom_index[chrom] = {"starts": starts, "entries": entries, "max_len": max_len}

    return gene_meta, chrom_index


def read_matches_strand(read: pysam.AlignedSegment, gene_strand: str) -> bool:
    # For STAR-aligned reads: read.is_reverse indicates alignment to reverse strand.
    # Gene strand '+' expects read NOT reverse; '-' expects read reverse.
    if gene_strand == "+":
        return not read.is_reverse
    return read.is_reverse


def genes_overlapping_psite(chrom_index, chrom: str, psite: int):
    """
    Return list of (gene_id, strand) where psite is within any CDS union interval.
    Uses a bounded backward scan from the last interval start <= psite.
    """
    if chrom not in chrom_index:
        return []

    starts = chrom_index[chrom]["starts"]   # list[int]
    entries = chrom_index[chrom]["entries"] # list[(s,e,gene,strand)]
    max_len = chrom_index[chrom]["max_len"]

    i = bisect_right(starts, int(psite)) - 1
    if i < 0:
        return []

    hits = []
    cutoff = int(psite) - int(max_len)

    j = i
    while j >= 0:
        s, e, gene, strand = entries[j]
        if s < cutoff:
            break
        if s <= psite < e:
            hits.append((gene, strand))
        j -= 1

    return hits


def is_frame0(psite: int, cds_start: int, strand: str) -> bool:
    """
    Frame relative to CDS start along the coding direction.
    '+' strand: distance = psite - cds_start
    '-' strand: distance = cds_start - psite
    """
    if strand == "+":
        return ((psite - cds_start) % 3) == 0
    return ((cds_start - psite) % 3) == 0


def main():
    args = parse_args()

    offsets = load_offsets(args.psite_offsets)
    gene_meta, chrom_index = parse_gtf_cds_union(args.gtf)

    gene_counts = defaultdict(int)
    stats = defaultdict(int)

    bam = pysam.AlignmentFile(args.bam, "rb")

    for read in bam.fetch(until_eof=True):
        stats["reads_seen"] += 1

        if read.is_unmapped:
            stats["unmapped"] += 1
            continue
        if read.is_secondary or read.is_supplementary:
            stats["secondary_or_supp"] += 1
            continue
        if (not args.include_duplicates) and read.is_duplicate:
            stats["duplicates_skipped"] += 1
            continue
        if read.mapping_quality < args.min_mapq:
            stats["low_mapq"] += 1
            continue

        psite, _ = get_psite(read, offsets)   # UPDATED: unpack (psite, read_len)
        if psite is None:
            stats["no_psite"] += 1
            continue

        chrom = bam.get_reference_name(read.reference_id)
        hits = genes_overlapping_psite(chrom_index, chrom, psite)
        if not hits:
            stats["psite_not_in_cds"] += 1
            continue

        # Strand-consistent hits only
        strand_hits = [(gene, strand) for gene, strand in hits if read_matches_strand(read, strand)]
        if not strand_hits:
            stats["strand_mismatch"] += 1
            continue

        # Ambiguity handling (multi-gene overlap)
        if args.skip_multi_overlap and len({g for g, _ in strand_hits}) > 1:
            stats["ambiguous_multi_gene"] += 1
            continue

        # Count
        counted_any = False
        for gene, strand in strand_hits:
            if args.inframe_only:
                cds_start = gene_meta[gene]["cds_start"]
                if not is_frame0(psite, cds_start, strand):
                    continue
            gene_counts[gene] += 1
            counted_any = True

        if counted_any:
            stats["counted"] += 1
        else:
            stats["filtered_inframe"] += 1

    bam.close()

    rows = []
    for g, meta in gene_meta.items():
        c = gene_counts.get(g, 0)
        cds_len = meta["cds_len"]
        gene_name = meta.get("gene_name", "")
        rows.append((args.sample, g, gene_name, c, cds_len))

    out_gene = pd.DataFrame(
        rows,
        columns=["sample", "gene_id", "gene_name", "cds_psite_count", "cds_length"],
    ).sort_values(["gene_id"])

    # TPM from CDS-length-normalized counts
    if not out_gene.empty:
        denom_kb = out_gene["cds_length"].replace(0, pd.NA) / 1000.0
        rpk = (out_gene["cds_psite_count"] / denom_kb).fillna(0.0)
        total_rpk = rpk.sum()
        out_gene["cds_tpm"] = (rpk / total_rpk * 1e6) if total_rpk > 0 else 0.0

    out_gene.to_csv(f"{args.out_prefix}.gene.tsv", sep="\t", index=False)

    pd.DataFrame(
        sorted(stats.items(), key=lambda x: x[0]),
        columns=["metric", "value"],
    ).to_csv(f"{args.out_prefix}.stats.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()