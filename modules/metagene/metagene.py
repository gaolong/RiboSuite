#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Metagene QC around start/stop using a P-site shifted BAM.")
    p.add_argument("--bam", required=True, help="P-site shifted BAM (read start represents P-site).")
    p.add_argument("--gtf", required=True, help="GTF annotation.")
    p.add_argument("--sample", required=True, help="Sample id prefix for outputs.")
    p.add_argument("--window", type=int, default=50, help="Window size (+/- N nt). Default: 50")
    p.add_argument("--max_transcripts_per_gene", type=int, default=1,
                   help="Keep up to N transcripts per gene (roughly: longest CDS). Default: 1")
    return p.parse_args()


def parse_gtf_attributes(attr_str: str) -> dict:
    # Standard GTF attributes: key "value";
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        d[m.group(1)] = m.group(2)
    return d


def load_cds_boundaries(gtf_path: str, max_tx_per_gene: int = 1):
    """
    Returns:
      start_sites[(chrom, strand)] -> set(genomic_position_0based)
      stop_sites[(chrom, strand)]  -> set(genomic_position_0based)

    We approximate:
      - "start site" = CDS start boundary (5' boundary of CDS)
      - "stop site"  = CDS end boundary (3' boundary of CDS, last CDS base)
    """
    # Collect CDS intervals per transcript
    tx_cds = defaultdict(list)  # (gene_id, tx_id, chrom, strand) -> list of (start0, end0_inclusive)
    tx_len = defaultdict(int)

    with open(gtf_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "CDS":
                continue

            a = parse_gtf_attributes(attrs)
            gene_id = a.get("gene_id")
            tx_id = a.get("transcript_id")
            if gene_id is None or tx_id is None:
                continue

            # GTF is 1-based inclusive; convert to 0-based inclusive
            s0 = int(start) - 1
            e0 = int(end) - 1

            key = (gene_id, tx_id, chrom, strand)
            tx_cds[key].append((s0, e0))
            tx_len[key] += (e0 - s0 + 1)

    # For each gene, pick up to N transcripts by total CDS length (fast heuristic)
    gene_to_txs = defaultdict(list)
    for (gene_id, tx_id, chrom, strand), cds_list in tx_cds.items():
        gene_to_txs[gene_id].append(((gene_id, tx_id, chrom, strand), tx_len[(gene_id, tx_id, chrom, strand)]))

    selected_keys = set()
    for gene_id, entries in gene_to_txs.items():
        entries.sort(key=lambda x: x[1], reverse=True)
        for (k, _l) in entries[:max_tx_per_gene]:
            selected_keys.add(k)

    start_sites = defaultdict(set)
    stop_sites = defaultdict(set)

    for key in selected_keys:
        gene_id, tx_id, chrom, strand = key
        cds_list = tx_cds[key]
        if not cds_list:
            continue

        # Determine CDS boundaries
        min_s = min(s for s, e in cds_list)
        max_e = max(e for s, e in cds_list)

        # Start site = 5' CDS boundary; Stop site = 3' CDS boundary
        if strand == "+":
            start_pos = min_s
            stop_pos = max_e
        else:
            # For minus strand, 5' is at higher coordinate; 3' at lower coordinate
            start_pos = max_e
            stop_pos = min_s

        start_sites[(chrom, strand)].add(start_pos)
        stop_sites[(chrom, strand)].add(stop_pos)

    return start_sites, stop_sites


def metagene_counts_from_bam(bam_path: str, start_sites, stop_sites, window: int):
    """
    BAM is assumed to be P-site shifted so that read.reference_start is the P-site (0-based).
    We aggregate counts by checking if psite - rel is in the site set for rel in [-window, +window].
    """
    n = 2 * window + 1
    start_counts = np.zeros(n, dtype=np.int64)
    stop_counts = np.zeros(n, dtype=np.int64)

    bam = pysam.AlignmentFile(bam_path, "rb")

    # Precompute rel positions array
    rels = np.arange(-window, window + 1, dtype=np.int64)

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        chrom = bam.get_reference_name(read.reference_id)
        strand = "-" if read.is_reverse else "+"
        key = (chrom, strand)

        psite = read.reference_start  # 0-based

        # START: if (psite - rel) is a start site, then this read contributes to rel
        sset = start_sites.get(key)
        if sset:
            # vectorized check: build candidate site positions = psite - rel
            candidates = (psite - rels).tolist()
            # loop over small window; window is small (e.g. 50), so this is fine
            for i, cand in enumerate(candidates):
                if cand in sset:
                    start_counts[i] += 1

        # STOP
        tset = stop_sites.get(key)
        if tset:
            candidates = (psite - rels).tolist()
            for i, cand in enumerate(candidates):
                if cand in tset:
                    stop_counts[i] += 1

    bam.close()
    return rels, start_counts, stop_counts


def write_tsv(sample: str, rels, counts, kind: str):
    df = pd.DataFrame({"rel_pos": rels, "count": counts})
    out = f"{sample}.metagene.{kind}.tsv"
    df.to_csv(out, sep="\t", index=False)
    return out, df


def plot_metagene(sample: str, start_df: pd.DataFrame, stop_df: pd.DataFrame, window: int):
    out_png = f"{sample}.metagene.png"

    plt.figure()
    plt.plot(start_df["rel_pos"], start_df["count"], label="Start (CDS boundary)")
    plt.plot(stop_df["rel_pos"], stop_df["count"], label="Stop (CDS boundary)")
    plt.axvline(0, linestyle="--")
    plt.xlabel("Position relative to site (nt)")
    plt.ylabel("P-site count")
    plt.title(f"Metagene QC (+/-{window} nt): {sample}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return out_png


def main():
    args = parse_args()

    start_sites, stop_sites = load_cds_boundaries(
        args.gtf,
        max_tx_per_gene=args.max_transcripts_per_gene
    )

    rels, start_counts, stop_counts = metagene_counts_from_bam(
        args.bam,
        start_sites,
        stop_sites,
        window=args.window
    )

    start_tsv, start_df = write_tsv(args.sample, rels, start_counts, "start")
    stop_tsv, stop_df = write_tsv(args.sample, rels, stop_counts, "stop")
    png = plot_metagene(args.sample, start_df, stop_df, args.window)

    print(f"[OK] Wrote: {start_tsv}, {stop_tsv}, {png}")


if __name__ == "__main__":
    main()
