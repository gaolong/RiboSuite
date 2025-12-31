#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt


# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Metagene QC with on-the-fly P-site calculation (sum-normalized)"
    )
    p.add_argument("--bam", required=True, help="Genome BAM (not P-site shifted)")
    p.add_argument("--gtf", required=True, help="GTF annotation")
    p.add_argument("--offset", required=True, help="P-site offset table (TSV)")
    p.add_argument("--sample", required=True, help="Sample ID prefix")
    p.add_argument("--window", type=int, default=50)
    p.add_argument("--min_cds_len", type=int, default=150)
    p.add_argument("--max_transcripts_per_gene", type=int, default=1)
    p.add_argument(
        "--require_canonical_codons",
        action="store_true",
        help="Require start_codon and stop_codon features",
    )
    return p.parse_args()


# --------------------------------------------------
# GTF parsing
# --------------------------------------------------
def parse_gtf_attributes(attr_str):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        d[m.group(1)] = m.group(2)
    return d


def load_cds_boundaries(gtf, min_cds_len, max_tx_per_gene, require_codons):
    tx_cds = defaultdict(list)
    tx_len = defaultdict(int)
    tx_has_start = defaultdict(bool)
    tx_has_stop = defaultdict(bool)

    with open(gtf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = line.rstrip().split("\t")
            a = parse_gtf_attributes(attrs)
            gene_id = a.get("gene_id")
            tx_id = a.get("transcript_id")
            if gene_id is None or tx_id is None:
                continue

            key = (gene_id, tx_id, chrom, strand)

            if feature == "CDS":
                s0 = int(start) - 1
                e0 = int(end) - 1
                tx_cds[key].append((s0, e0))
                tx_len[key] += (e0 - s0 + 1)
            elif feature == "start_codon":
                tx_has_start[key] = True
            elif feature == "stop_codon":
                tx_has_stop[key] = True

    gene_to_txs = defaultdict(list)
    for k, cds_len in tx_len.items():
        if cds_len < min_cds_len:
            continue
        if require_codons and not (tx_has_start[k] and tx_has_stop[k]):
            continue
        gene_to_txs[k[0]].append((k, cds_len))

    selected = set()
    for gene, entries in gene_to_txs.items():
        entries.sort(key=lambda x: x[1], reverse=True)
        for k, _ in entries[:max_tx_per_gene]:
            selected.add(k)

    start_sites = defaultdict(set)
    stop_sites = defaultdict(set)

    for gene_id, tx_id, chrom, strand in selected:
        cds_list = tx_cds[(gene_id, tx_id, chrom, strand)]
        min_s = min(s for s, _ in cds_list)
        max_e = max(e for _, e in cds_list)

        if strand == "+":
            start_pos = min_s
            stop_pos = max_e
        else:
            start_pos = max_e
            stop_pos = min_s

        start_sites[(chrom, strand)].add(start_pos)
        stop_sites[(chrom, strand)].add(stop_pos)

    return start_sites, stop_sites


# --------------------------------------------------
# Offset table
# --------------------------------------------------
def load_offsets(path):
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["read_length"], df["psite_offset"]))


# --------------------------------------------------
# Metagene core (COUNTS)
# --------------------------------------------------
def metagene_counts_from_bam(bam_path, start_sites, stop_sites, offset_table, window):
    n = 2 * window + 1
    rels = np.arange(-window, window + 1, dtype=np.int64)

    start_counts = np.zeros(n, dtype=np.int64)
    stop_counts = np.zeros(n, dtype=np.int64)

    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        L = read.query_length
        offset = offset_table.get(L)
        if offset is None:
            continue

        chrom = bam.get_reference_name(read.reference_id)
        strand = "-" if read.is_reverse else "+"
        key = (chrom, strand)

        if read.is_reverse:
            psite = read.reference_end - 1 - offset
        else:
            psite = read.reference_start + offset

        if strand == "+":
            candidates = psite - rels
        else:
            candidates = psite + rels

        sset = start_sites.get(key)
        if sset:
            for i, c in enumerate(candidates):
                if c in sset:
                    start_counts[i] += 1

        tset = stop_sites.get(key)
        if tset:
            for i, c in enumerate(candidates):
                if c in tset:
                    stop_counts[i] += 1

    bam.close()
    return rels, start_counts, stop_counts


# --------------------------------------------------
# Output
# --------------------------------------------------
def write_tsv(sample, rels, counts, kind):
    df = pd.DataFrame({"rel_pos": rels, "count": counts})
    out = f"{sample}.metagene.{kind}.tsv"
    df.to_csv(out, sep="\t", index=False)
    return df


def plot_metagene(sample, rels, start_counts, stop_counts, window):
    out_png = f"{sample}.metagene.png"

    # global sum normalization (ribowaltz-like frequency)
    start_freq = start_counts / start_counts.sum() if start_counts.sum() > 0 else start_counts
    stop_freq = stop_counts / stop_counts.sum() if stop_counts.sum() > 0 else stop_counts

    plt.figure(figsize=(6, 4))
    plt.plot(rels, start_freq, label="Start")
    plt.plot(rels, stop_freq, label="Stop")
    plt.axvline(0, linestyle="--", color="black", linewidth=0.8)
    plt.xlabel("Position relative to site (nt)")
    plt.ylabel("P-site frequency")
    plt.title(f"Metagene QC (+/-{window} nt): {sample}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return out_png


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    args = parse_args()

    start_sites, stop_sites = load_cds_boundaries(
        args.gtf,
        args.min_cds_len,
        args.max_transcripts_per_gene,
        args.require_canonical_codons,
    )

    offset_table = load_offsets(args.offset)

    rels, start_counts, stop_counts = metagene_counts_from_bam(
        args.bam,
        start_sites,
        stop_sites,
        offset_table,
        args.window,
    )

    write_tsv(args.sample, rels, start_counts, "start")
    write_tsv(args.sample, rels, stop_counts, "stop")

    png = plot_metagene(
        args.sample,
        rels,
        start_counts,
        stop_counts,
        args.window,
    )

    print(f"[OK] Metagene QC completed: {png}")


if __name__ == "__main__":
    main()
