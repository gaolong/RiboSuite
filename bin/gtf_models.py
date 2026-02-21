#!/usr/bin/env python3
"""
CDS-only model loader.

- Does not parse start_codon/stop_codon features.
- Defines ORF start/stop positions using CDS boundaries.
- Derives UTR intervals in tx coords as: exon_tx - cds_tx.

Conventions:
- Genomic and tx coordinates are 0-based half-open [start, end)
- anchor_tx = cds_start_tx (no phase adjustment)
"""

from collections import defaultdict
from typing import Dict, List, Tuple, Any, Optional

from psite_utils import parse_attrs, genomic_to_tx_spliced
from intervals import merge_intervals, subtract_intervals


def keep_longest_tx_per_gene(models: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    best = {}
    for m in models:
        gene = m["gene"]
        cds_len = m.get("cds_len", 0)
        if gene not in best or cds_len > best[gene]["cds_len"]:
            best[gene] = {"model": m, "cds_len": cds_len}
    return [v["model"] for v in best.values()]


def _build_exons_txpos(exons: List[Tuple[int, int]]) -> Tuple[List[Tuple[int, int, int]], int]:
    exons.sort(key=lambda x: x[0])
    exons_txpos: List[Tuple[int, int, int]] = []
    cursor = 0
    for s, e in exons:
        exons_txpos.append((s, e, cursor))
        cursor += (e - s)
    return exons_txpos, cursor


def _tx_intervals_from_genomic_blocks(
    blocks: List[Tuple[int, int]],
    exons_txpos: List[Tuple[int, int, int]],
) -> List[Tuple[int, int]]:
    """
    Convert genomic blocks [s,e) into tx intervals [tx_s, tx_e), potentially split by splicing.

    Works even if a block spans introns: it will map only the exon-overlapping parts.
    """
    out: List[Tuple[int, int]] = []
    for bs, be in blocks:
        # intersect with each exon; map overlap to tx
        for es, ee, tx_s0 in exons_txpos:
            os = max(bs, es)
            oe = min(be, ee)
            if os < oe:
                tx_s = tx_s0 + (os - es)
                tx_e = tx_s0 + (oe - es)
                out.append((tx_s, tx_e))
    return merge_intervals(out)


def load_tx_models_cds_only(
    gtf: str,
    min_cds_len: int = 0,
    protein_coding_only: bool = False,
) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """
    Build models using exon + CDS features only.

    protein_coding_only:
      - True for GENCODE-like transcript QC
      - False for ORF GTFs (recommended)
    """

    tx_exons = defaultdict(list)  # key -> [(s,e)]
    tx_cds = defaultdict(list)    # key -> [(s,e)]
    gene_type = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue

            chrom, _src, feature, start, end, _score, strand, _phase, attr = fields
            attrs = parse_attrs(attr)

            if feature == "gene":
                gid = attrs.get("gene_id")
                if gid:
                    gene_type[gid] = attrs.get("gene_type") or attrs.get("gene_biotype")
                continue

            gene = attrs.get("gene_id")
            tx = attrs.get("transcript_id")
            if not gene or not tx:
                continue

            if protein_coding_only and gene_type.get(gene) != "protein_coding":
                continue

            start0 = int(start) - 1
            end0 = int(end)
            key = (gene, tx, chrom, strand)

            if feature == "exon":
                tx_exons[key].append((start0, end0))
            elif feature == "CDS":
                tx_cds[key].append((start0, end0))

    models: List[Dict[str, Any]] = []
    tx_meta: Dict[str, Dict[str, Any]] = {}

    for key, exons in tx_exons.items():
        gene, tx, chrom, strand = key
        cds_blocks = tx_cds.get(key)
        if not cds_blocks:
            continue

        cds_len = sum(e - s for s, e in cds_blocks)
        if cds_len < min_cds_len:
            continue

        exons_txpos, tx_len = _build_exons_txpos(exons)

        # Convert exon span into tx intervals: by definition it is [0, tx_len)
        exon_tx_intervals = [(0, tx_len)]

        # Convert CDS genomic blocks -> tx intervals (spliced)
        cds_tx_intervals = _tx_intervals_from_genomic_blocks(cds_blocks, exons_txpos)
        if not cds_tx_intervals:
            # CDS blocks didn't overlap exons in a mappable way
            continue

        # Translation bounds in tx coords = union CDS span
        cds_lo_tx = min(s for s, _ in cds_tx_intervals)
        cds_hi_tx = max(e for _, e in cds_tx_intervals)

        # Anchor frame at CDS start (no phase)
        # For minus strand, this still anchors at the "leftmost in tx coords" CDS start,
        # which is consistent with using tx coords and frame formula that accounts for strand.
        anchor_tx = cds_lo_tx

        # Derive UTR intervals = exon - CDS (all in tx coords)
        utr_tx_intervals = subtract_intervals(exon_tx_intervals, cds_tx_intervals)

        fetch_start = min(s for s, e in exons)
        fetch_end = max(e for s, e in exons)

        model = {
            "gene": gene,
            "tx": tx,
            "chrom": chrom,
            "strand": strand,
            "exons": exons_txpos,

            "cds_blocks_g": cds_blocks,         # genomic
            "cds_tx": cds_tx_intervals,         # tx coords, spliced
            "utr_tx": utr_tx_intervals,         # derived, tx coords

            "cds_lo_tx": cds_lo_tx,
            "cds_hi_tx": cds_hi_tx,
            "anchor_tx": anchor_tx,

            "fetch_start": fetch_start,
            "fetch_end": fetch_end,
            "cds_len": cds_len,
            "tx_len": tx_len,
        }
        models.append(model)

        tx_meta[tx] = {
            "gene": gene,
            "chrom": chrom,
            "strand": strand,
            "cds_tx": cds_tx_intervals,
            "utr_tx": utr_tx_intervals,
            "cds_lo_tx": cds_lo_tx,
            "cds_hi_tx": cds_hi_tx,
            "anchor_tx": anchor_tx,
            "cds_len": cds_len,
            "tx_len": tx_len,
            "fetch_start": fetch_start,
            "fetch_end": fetch_end,
        }

    return models, tx_meta