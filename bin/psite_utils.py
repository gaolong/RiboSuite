#!/usr/bin/env python3
"""
psite_utils.py

Shared helpers for RiboSuite modules.

Conventions
-----------
- Genomic and transcript coordinates are 0-based.
- Intervals are half-open: [start, end)
- P-site returned is a genomic 0-based coordinate.

Design notes
------------
- get_psite() can use read orientation (read.is_reverse) OR be forced to use a
  provided strand ('+' / '-') to keep consistency with transcript models.
"""

from __future__ import annotations

import re
from typing import Dict, Optional, Tuple, List, Iterable

import pandas as pd


# --------------------------------------------------
# offsets
# --------------------------------------------------
def load_offsets(path: str) -> Dict[int, int]:
    """
    Load read_length -> psite_offset mapping from TSV.

    Accepts either column name: psite_offset or offset.

    Returns
    -------
    dict[int, int]
        Keys: read_length (int)
        Values: psite_offset (int)
    """
    df = pd.read_csv(path, sep="\t")
    col = "psite_offset" if "psite_offset" in df.columns else "offset"
    return dict(zip(df["read_length"].astype(int), df[col].astype(int)))


# --------------------------------------------------
# P-site
# --------------------------------------------------
def get_psite(
    read,
    offsets: Dict[int, int],
    strand: Optional[str] = None,
) -> Tuple[Optional[int], Optional[int]]:
    """
    Infer P-site genomic position (0-based) using length-specific offsets.

    Parameters
    ----------
    read : pysam.AlignedSegment
    offsets : dict[int, int]
        read_length -> psite_offset
    strand : Optional[str]
        If provided, must be '+' or '-'. This overrides read.is_reverse and forces
        P-site placement to follow the feature/transcript strand.

    Returns
    -------
    (psite, read_len) or (None, None)
    """
    read_len = read.query_length
    off = offsets.get(read_len)
    if off is None:
        return None, None

    if strand is None:
        is_reverse = read.is_reverse
    else:
        if strand not in {"+", "-"}:
            raise ValueError(f"strand must be '+' or '-', got {strand!r}")
        is_reverse = (strand == "-")

    # reference_end can be None in rare edge cases; be defensive
    if is_reverse:
        if read.reference_end is None:
            return None, None
        psite = read.reference_end - off - 1
    else:
        psite = read.reference_start + off

    return psite, read_len


# --------------------------------------------------
# intervals
# --------------------------------------------------
def in_interval(pos: int, intervals: Iterable[Tuple[int, int]]) -> bool:
    """
    Return True if pos lies in any interval in a list of half-open intervals [s,e).
    """
    for s, e in intervals:
        if s <= pos < e:
            return True
    return False


# --------------------------------------------------
# GTF helpers
# --------------------------------------------------
_ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attrs(attr: str) -> Dict[str, str]:
    """
    Parse a GTF attributes field into a dict.

    Example:
      gene_id "ENSG..."; transcript_id "ENST..."; gene_type "protein_coding";
    """
    return dict(_ATTR_RE.findall(attr))


def genomic_to_tx_spliced(
    pos: int,
    exons_txpos: List[Tuple[int, int, int]],
) -> Optional[int]:
    """
    Map a genomic 0-based position to spliced transcript coordinate (0-based).

    Parameters
    ----------
    pos : int
        Genomic position (0-based).
    exons_txpos : list of (g_start, g_end, tx_start)
        Exons in genomic order. Each exon maps:
          [g_start, g_end) -> [tx_start, tx_start + (g_end - g_start))

    Returns
    -------
    tx_pos : int or None
        Transcript coordinate if pos falls within an exon; otherwise None.
    """
    for g_s, g_e, tx_s in exons_txpos:
        if g_s <= pos < g_e:
            return tx_s + (pos - g_s)
    return None


