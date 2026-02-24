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
We compute an anchor-relative coordinate (rel) in transcript (spliced) space:
  + strand: rel = tx_pos - anchor_tx
  - strand: rel = anchor_tx - tx_pos
Frame is then rel % 3 and codon index is rel // 3.

Update (Uniformity)
-------------------
Uniformity is computed as the fraction of CDS codons whose per-codon periodicity
in translating position (frame 0 relative to anchor) exceeds a threshold.

Update (Drop-off)
-----------------
We compute a stop-codon drop-off metric using a symmetric window of N codons
around the CDS end, counting only translating-frame (frame 0) P-sites:

  upstream   = frame0 P-sites in last N CDS codons
  downstream = frame0 P-sites in first N codons immediately after CDS end
  dropoff_raw = upstream / (upstream + downstream)

Because dropoff_raw collapses low-support extremes (e.g., up=1,down=0 vs up=100,down=0),
we also report:
  dropoff_den_f0 = upstream + downstream   (evidence)
  dropoff_smoothed = (upstream + alpha) / (upstream + downstream + alpha + beta)

We optionally suppress dropoff_raw when evidence is too low (denom < min_dropoff_den_f0),
setting dropoff_status="low_support".
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

DEFAULT_DROPOFF_WINDOW_CODONS = 10
DEFAULT_DROPOFF_ALPHA = 1.0
DEFAULT_DROPOFF_BETA = 10.0
DEFAULT_MIN_DROPOFF_DEN_F0 = 5


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
    p.add_argument(
        "--min_len_codons",
        type=int,
        default=30,
        help="Minimum CDS length in codons (applied as min_cds_len = codons*3)",
    )
    p.add_argument("--min_frame_fraction", type=float, default=0.7)
    p.add_argument("--mapq", type=int, default=100)

    # Uniformity metric settings (metric computed regardless of rule usage)
    p.add_argument(
        "--uniformity_threshold",
        type=float,
        default=DEFAULT_UNIFORMITY_THRESHOLD,
        help="Per-codon periodicity threshold for uniformity (default: 1/3)",
    )

    # Drop-off metric settings
    p.add_argument(
        "--dropoff_window_codons",
        type=int,
        default=DEFAULT_DROPOFF_WINDOW_CODONS,
        help="Window size (codons) for drop-off around CDS end (default: 10). "
             "Counts frame0 P-sites in last N CDS codons vs first N codons after CDS end.",
    )
    p.add_argument(
        "--dropoff_alpha",
        type=float,
        default=DEFAULT_DROPOFF_ALPHA,
        help="Pseudo-count alpha for dropoff smoothing (default: 1.0).",
    )
    p.add_argument(
        "--dropoff_beta",
        type=float,
        default=DEFAULT_DROPOFF_BETA,
        help="Pseudo-count beta for dropoff smoothing (default: 1.0).",
    )
    p.add_argument(
        "--min_dropoff_den_f0",
        type=int,
        default=DEFAULT_MIN_DROPOFF_DEN_F0,
        help="Minimum (up+down) frame0 P-sites required to report dropoff_raw (default: 5). "
             "dropoff_smoothed is still reported when denom>0.",
    )

    # Recommended default for ORF GTFs is False; for pure GENCODE QC set True.
    p.add_argument(
        "--protein_coding_only",
        action="store_true",
        help="Filter to protein_coding genes based on gene_type/gene_biotype in GTF",
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
def assign_psites_to_model(
    bam,
    offsets,
    m,
    mapq,
    uniformity_threshold=DEFAULT_UNIFORMITY_THRESHOLD,
    dropoff_window_codons=DEFAULT_DROPOFF_WINDOW_CODONS,
    dropoff_alpha=DEFAULT_DROPOFF_ALPHA,
    dropoff_beta=DEFAULT_DROPOFF_BETA,
    min_dropoff_den_f0=DEFAULT_MIN_DROPOFF_DEN_F0,
):
    """
    Count P-sites mapped into CDS, by frame (tx_pos anchored at anchor_tx),
    plus codon coverage, uniformity in CDS space, and stop-codon drop-off.

    Frame definition (consistent with frame_periodicity):
      + strand: frame = (tx_pos - anchor_tx) % 3
      - strand: frame = (anchor_tx - tx_pos) % 3

    Drop-off bins are defined in translation-direction codon coordinates:
      codon = rel // 3  where rel is anchor-relative (>=0).

    Metrics
    -------
    - total/fcounts: CDS-only P-sites by frame
    - coverage_fraction: fraction of CDS codons with >=1 P-site (any frame)
    - uniformity: fraction of CDS codons where (frame0 / all_frames_in_codon) > threshold
    - dropoff_up_f0: frame0 P-sites in last N CDS codons
    - dropoff_down_f0: frame0 P-sites in first N codons after CDS end
    - dropoff_den_f0: evidence (up+down)
    - dropoff_raw: up/(up+down) if evidence >= min_dropoff_den_f0 else None
    - dropoff_smoothed: (up+alpha)/(up+down+alpha+beta) if denom>0 else None
    - dropoff_status: ok/no_signal/no_upstream/no_downstream/low_support
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

    # codon index (anchor-relative) -> [pos0,pos1,pos2] counts (CDS only)
    codon_pos_counts = defaultdict(lambda: [0, 0, 0])

    # Drop-off counters (frame0 only)
    up_f0 = 0
    down_f0 = 0

    # Define drop-off codon windows (anchor-relative codon indices)
    N = max(0, int(dropoff_window_codons))
    up_start = max(0, cds_len_codons - N)
    up_end = cds_len_codons
    down_start = cds_len_codons
    down_end = cds_len_codons + N

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

        # Anchor-relative coordinate (translation direction)
        if strand == "+":
            rel = tx_pos - anchor_tx
        else:
            rel = anchor_tx - tx_pos

        # Exclude upstream of start (5'UTR direction)
        if rel < 0:
            continue

        frame = rel % 3
        codon = rel // 3

        # -----------------------
        # Drop-off (frame0 only)
        # -----------------------
        if N > 0 and frame == 0:
            if up_start <= codon < up_end:
                up_f0 += 1
            elif down_start <= codon < down_end:
                down_f0 += 1

        # -----------------------
        # CDS-only metrics below
        # -----------------------
        cds_pos = tx_pos_to_cds_pos(tx_pos, cds_tx)
        if cds_pos is None:
            continue

        fcounts[frame] += 1
        total += 1

        if codon < cds_len_codons:
            codon_pos_counts[codon][frame] += 1

    # Coverage fraction
    codons_with_psite = sum(1 for v in codon_pos_counts.values() if sum(v) > 0)
    coverage_fraction = (codons_with_psite / cds_len_codons) if cds_len_codons > 0 else 0.0

    # Uniformity
    good = 0
    for v in codon_pos_counts.values():
        tot = v[0] + v[1] + v[2]
        if tot == 0:
            continue
        if (v[0] / tot) > uniformity_threshold:
            good += 1

    uniformity = (good / cds_len_codons) if cds_len_codons > 0 else 0.0

    # Drop-off
    denom = up_f0 + down_f0

    if denom == 0:
        dropoff_status = "no_signal"
    elif up_f0 == 0:
        dropoff_status = "no_upstream"
    elif down_f0 == 0:
        dropoff_status = "no_downstream"
    else:
        dropoff_status = "ok"

    # Raw ratio only if enough evidence
    if denom >= min_dropoff_den_f0:
        dropoff_raw = up_f0 / denom
    else:
        dropoff_raw = None
        if denom > 0:
            dropoff_status = "low_support"

    # Smoothed ratio (helpful even when downstream==0, and distinguishes low vs high support)
    dropoff_smoothed = (
        (up_f0 + dropoff_alpha) / (denom + dropoff_alpha + dropoff_beta)
        if denom > 0 else None
    )

    return (
        total,
        fcounts,
        coverage_fraction,
        uniformity,
        dropoff_raw,
        dropoff_smoothed,
        up_f0,
        down_f0,
        denom,
        dropoff_status,
    )


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

    models = load_tx_models(
        args.gtf,
        min_cds_len=args.min_len_codons * 3,
        protein_coding_only=args.protein_coding_only,
    )
    if not models:
        raise RuntimeError("No transcript models loaded (check GTF features/attrs)")

    rows = []

    for m in models:
        (
            total,
            fcounts,
            cov_frac,
            unif,
            dropoff_raw,
            dropoff_smoothed,
            up_f0,
            down_f0,
            drop_den,
            drop_status,
        ) = assign_psites_to_model(
            bam,
            offsets,
            m,
            args.mapq,
            uniformity_threshold=args.uniformity_threshold,
            dropoff_window_codons=args.dropoff_window_codons,
            dropoff_alpha=args.dropoff_alpha,
            dropoff_beta=args.dropoff_beta,
            min_dropoff_den_f0=args.min_dropoff_den_f0,
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

            # Drop-off schema (ratio + evidence + status)
            "dropoff_window_codons": args.dropoff_window_codons,
            "dropoff_up_f0": up_f0,
            "dropoff_down_f0": down_f0,
            "dropoff_den_f0": drop_den,
            "dropoff_status": drop_status,
            "dropoff_raw": (round(dropoff_raw, 3) if dropoff_raw is not None else None),
            "dropoff_smoothed": (round(dropoff_smoothed, 3) if dropoff_smoothed is not None else None),

            "rules_passed": ",".join(rules_passed),
            "rules_failed": ",".join(rules_failed),
        })

    df = pd.DataFrame(rows)
    out_tsv = f"{args.out_prefix}.translated_orfs.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"[translated_orfs] written {out_tsv}")


if __name__ == "__main__":
    main()