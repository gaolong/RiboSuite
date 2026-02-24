#!/usr/bin/env python3
"""
call_translated_orfs.py

Rule-based identification of translated ORFs using CDS-only frame periodicity
(no start/stop-codon dependency).

This version reuses shared utilities in bin/:
- bin.gtf_models.load_tx_models
- bin.psite_utils.load_offsets, get_psite, genomic_to_tx_spliced

Key idea
--------
We compute a CDS-only coordinate (cds_pos) by accumulating across the CDS
intervals in transcript coordinates. Frame is then cds_pos % 3.

Update (Uniformity)
-------------------
In addition to codon coverage ("coverage_fraction"), we compute "uniformity":

  For each codon, compute the fraction of P-sites that land in the translating
  position (frame 0 relative to anchor) among all P-sites in that codon.
  Count codons where this fraction > threshold (default 1/3), and divide by
  total codons in ORF.

  uniformity = (# codons with per-codon periodicity > threshold) / (total codons)

NOTE: In this implementation, the translating position is defined as frame 0
relative to anchor_tx (i.e., rel % 3 == 0).
"""

import argparse
from collections import defaultdict

import pysam
import pandas as pd

# --------------------------------------------------
# Ensure imports from project root work
# --------------------------------------------------
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../../"))
sys.path.insert(0, ROOT_DIR)

from bin.gtf_models import load_tx_models
from bin.psite_utils import load_offsets, get_psite, genomic_to_tx_spliced


# --------------------------------------------------
# constants
# --------------------------------------------------
MIN_COVERAGE_FRACTION = 0.25
DEFAULT_UNIFORMITY_THRESHOLD = 1.0 / 3.0


# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Rule-based identification of translated ORFs "
                    "using CDS-only frame periodicity (no start/stop-codon dependency)"
    )
    p.add_argument("--bam", required=True, help="Deduplicated, sorted BAM")
    p.add_argument("--gtf", required=True, help="GTF annotation")
    p.add_argument("--psite_offsets", required=True, help="P-site offset table (TSV)")
    p.add_argument("--out_prefix", required=True)

    p.add_argument("--min_psites", type=int, default=30)
    p.add_argument("--min_len_codons", type=int, default=30,
                   help="Minimum CDS length in codons (applied as min_cds_len = codons*3)")
    p.add_argument("--min_frame_fraction", type=float, default=0.7)
    p.add_argument("--mapq", type=int, default=100)

    # Uniformity metric settings (metric computed regardless of rule usage)
    p.add_argument(
        "--uniformity_threshold",
        type=float,
        default=DEFAULT_UNIFORMITY_THRESHOLD,
        help="Per-codon periodicity threshold for uniformity (default: 1/3)"
    )

    # Recommended default for ORF GTFs is False; for pure GENCODE QC set True.
    p.add_argument(
        "--protein_coding_only",
        action="store_true",
        help="Filter to protein_coding genes based on gene_type/gene_biotype in GTF"
    )

    return p.parse_args()


# --------------------------------------------------
# helpers: tx_pos -> CDS-only coordinate
# --------------------------------------------------
def cds_len_from_intervals(cds_tx_intervals):
    """Total CDS length in nt from CDS tx intervals."""
    return sum(e - s for s, e in cds_tx_intervals)


def tx_pos_to_cds_pos(tx_pos, cds_tx_intervals):
    """
    Convert transcript coordinate -> CDS-only coordinate by accumulating lengths
    across CDS tx intervals (half-open [s,e)).

    Returns:
      cds_pos (0-based) if tx_pos is in CDS; otherwise None
    """
    acc = 0
    for s, e in cds_tx_intervals:
        if tx_pos < s:
            break
        if s <= tx_pos < e:
            return acc + (tx_pos - s)
        acc += (e - s)
    return None


# --------------------------------------------------
# core logic
# --------------------------------------------------
def assign_psites_to_model(bam, offsets, m, mapq, uniformity_threshold=DEFAULT_UNIFORMITY_THRESHOLD):
    """
    Count P-sites mapped into CDS, by frame (tx_pos anchored at anchor_tx),
    plus codon coverage and uniformity in CDS space.

    Frame definition (consistent with frame_periodicity):
      + strand: frame = (tx_pos - anchor_tx) % 3
      - strand: frame = (anchor_tx - tx_pos) % 3

    Notes
    -----
    - Requires m["anchor_tx"] to represent the translation start position in tx coords
      (strand-aware start codon proxy). Fix in gtf_models:
        anchor_tx = cds_lo_tx if strand == "+" else (cds_hi_tx - 1)
    - Coverage is computed as fraction of CDS codons with >=1 P-site (any frame).
    - Uniformity is computed by per-codon 3-nt periodicity:
        periodicity_codon = (P-sites in translating position) / (P-sites in codon)
      where translating position is frame 0 relative to anchor.
    """
    fcounts = [0, 0, 0]
    total = 0

    chrom = m["chrom"]
    strand = m["strand"]
    anchor_tx = m["anchor_tx"]

    max_offset = max(offsets.values())
    fs = max(0, m["fetch_start"] - max_offset)
    fe = m["fetch_end"] + max_offset

    cds_tx = m["cds_tx"]  # list of (s,e) in tx coords (spliced)
    cds_len_nt = cds_len_from_intervals(cds_tx)
    cds_len_codons = cds_len_nt // 3

    # codon index (anchor-relative) -> [pos0,pos1,pos2] counts
    codon_pos_counts = defaultdict(lambda: [0, 0, 0])

    for read in bam.fetch(chrom, fs, fe):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < mapq:
            continue

        psite, _rl = get_psite(read, offsets, strand=strand)
        if psite is None:
            continue

        tx_pos = genomic_to_tx_spliced(psite, m["exons"])
        if tx_pos is None:
            continue

        # Must be in CDS for translated-ORF rules
        cds_pos = tx_pos_to_cds_pos(tx_pos, cds_tx)
        if cds_pos is None:
            continue

        # Anchor-relative coordinate (periodicity axis)
        if strand == "+":
            rel = tx_pos - anchor_tx
        else:
            rel = anchor_tx - tx_pos

        # If anchor_tx is correct CDS start proxy, rel should be >= 0 for CDS reads.
        if rel < 0:
            continue

        frame = rel % 3
        codon = rel // 3

        fcounts[frame] += 1
        total += 1

        # Only tally codons within CDS length
        if codon < cds_len_codons:
            codon_pos_counts[codon][frame] += 1

    # Coverage fraction: codons with any P-site (any position)
    codons_with_psite = sum(1 for v in codon_pos_counts.values() if sum(v) > 0)
    coverage_fraction = (
        codons_with_psite / cds_len_codons
        if cds_len_codons > 0 else 0.0
    )

    # Uniformity: codons where per-codon periodicity in translating position > threshold
    # Translating position is frame 0 relative to anchor.
    good = 0
    for v in codon_pos_counts.values():
        tot = v[0] + v[1] + v[2]
        if tot == 0:
            continue
        periodicity = v[0] / tot
        if periodicity > uniformity_threshold:
            good += 1

    uniformity = (
        good / cds_len_codons
        if cds_len_codons > 0 else 0.0
    )

    return total, fcounts, coverage_fraction, uniformity


# --------------------------------------------------
# rules
# --------------------------------------------------
def rule_min_psites(total, min_psites):
    return total >= min_psites


def rule_frame_dominance(fcounts, min_fraction):
    s = sum(fcounts)
    if s == 0:
        return False, 0.0
    frac = max(fcounts) / s
    return frac >= min_fraction, frac


def rule_coverage_fraction(cov_frac, min_cov):
    return cov_frac >= min_cov


# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()

    offsets = load_offsets(args.psite_offsets)
    if not offsets:
        raise RuntimeError("No offsets loaded")

    bam = pysam.AlignmentFile(args.bam, "rb")

    # Reuse shared transcript/CDS models
    models = load_tx_models(
        args.gtf,
        min_cds_len=args.min_len_codons * 3,
        protein_coding_only=args.protein_coding_only,
    )

    if not models:
        raise RuntimeError("No transcript models loaded (check GTF features/attrs)")

    rows = []

    for m in models:
        total, fcounts, cov_frac, unif = assign_psites_to_model(
            bam, offsets, m, args.mapq, uniformity_threshold=args.uniformity_threshold
        )

        # hard output filter
        if total < args.min_psites:
            continue

        rules_passed = []
        rules_failed = []

        if rule_min_psites(total, args.min_psites):
            rules_passed.append("min_psites")
        else:
            rules_failed.append("min_psites")

        ok_frame, frac = rule_frame_dominance(fcounts, args.min_frame_fraction)
        if ok_frame:
            rules_passed.append("frame_dominance")
        else:
            rules_failed.append("frame_dominance")

        if rule_coverage_fraction(cov_frac, MIN_COVERAGE_FRACTION):
            rules_passed.append("coverage_fraction")
        else:
            rules_failed.append("coverage_fraction")

        rows.append({
            "gene_id": m["gene"],
            "gene_name": m.get("gene_name"),
            "transcript_id": m["tx"],
            "psite_total": total,
            "psite_f0": fcounts[0],
            "psite_f1": fcounts[1],
            "psite_f2": fcounts[2],
            "frame_fraction": round(frac, 3),
            "coverage_fraction": round(cov_frac, 3),
            "uniformity": round(unif, 3),
            "rules_passed": ",".join(rules_passed),
            "rules_failed": ",".join(rules_failed),
        })

    df = pd.DataFrame(rows)
    out_tsv = f"{args.out_prefix}.translated_orfs.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"[translated_orfs] written {out_tsv}")


if __name__ == "__main__":
    main()