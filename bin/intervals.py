#!/usr/bin/env python3
"""
Interval utilities for half-open intervals [start, end).
"""

from typing import List, Tuple


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping/adjacent half-open intervals.
    """
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:  # overlap or adjacency
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def subtract_intervals(
    a: List[Tuple[int, int]],
    b: List[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """
    Return A \ B for half-open intervals, assuming nothing about sorting.

    Example:
      A = [(0, 100)]
      B = [(10, 20), (30, 40)]
      -> [(0, 10), (20, 30), (40, 100)]
    """
    a = merge_intervals(a)
    b = merge_intervals(b)

    out: List[Tuple[int, int]] = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            bs, be = b[k]
            if bs > cur:
                out.append((cur, min(bs, e)))
            cur = max(cur, be)
            if cur >= e:
                break
            k += 1
        if cur < e:
            out.append((cur, e))
    return out