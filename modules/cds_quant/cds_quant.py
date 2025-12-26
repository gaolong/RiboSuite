#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
import re
from collections import defaultdict
from psite_utils import get_psite


def parse_args():
    p = argparse.ArgumentParser(
        description="CDS-level P-site quantification (on-the-fly)"
    )
    p.add_argument("--bam", required=True, help="STAR-aligned BAM")
    p.add_argument("--gtf", required=True, help="GTF annotation")
    p.add_argument("--offsets", required=True, help="P-site offsets TSV")
    p.add_argument("--sample", required=True, help="Sample ID")
    p.add_argument("--inframe_only", action="store_true",
                   help="Only count frame-0 P-sites")
    p.add_argument("--out_prefix", default="cds_quant")
    return p.parse_args()


def load_offsets(offset_file):
    """
    TSV with columns: read_length <tab> offset
    """
    df = pd.read_csv(offset_file, sep="\t", header=None)
    return dict(zip(df[0], df[1]))


def parse_gtf(gtf_file):
    cds_by_tx = defaultdict(list)
    gene_by_tx = {}

    attr_re = re.compile(r'(\S+) "([^"]+)"')

    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            chrom, _, feature, start, end, _, strand, frame, attr = line.strip().split("\t")
            if feature != "CDS":
                continue

            attrs = dict(attr_re.findall(attr))
            tx = attrs.get("transcript_id")
            gene = attrs.get("gene_id")

            if not tx or not gene:
                continue

            start = int(start) - 1
            end = int(end)

            cds_by_tx[tx].append((chrom, start, end, strand, frame))
            gene_by_tx[tx] = gene

    return cds_by_tx, gene_by_tx


def main():
    args = parse_args()

    offsets = load_offsets(args.offsets)
    cds_by_tx, gene_by_tx = parse_gtf(args.gtf)

    bam = pysam.AlignmentFile(args.bam, "rb")

    tx_counts = defaultdict(int)
    gene_counts = defaultdict(int)

    for tx, cds_list in cds_by_tx.items():
        for chrom, start, end, strand, frame in cds_list:
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped:
                    continue

                psite = get_psite(read, offsets)
                if psite is None or not (start <= psite < end):
                    continue

                if strand == "+" and read.is_reverse:
                    continue
                if strand == "-" and not read.is_reverse:
                    continue

                if args.inframe_only:
                    cds_frame = int(frame) if frame.isdigit() else 0
                    if (psite - start - cds_frame) % 3 != 0:
                        continue

                tx_counts[tx] += 1
                gene_counts[gene_by_tx[tx]] += 1

    bam.close()

    pd.DataFrame(
        [(args.sample, tx, c) for tx, c in tx_counts.items()],
        columns=["sample", "transcript_id", "cds_psite_count"]
    ).to_csv(f"{args.out_prefix}.transcript.tsv", sep="\t", index=False)

    pd.DataFrame(
        [(args.sample, g, c) for g, c in gene_counts.items()],
        columns=["sample", "gene_id", "cds_psite_count"]
    ).to_csv(f"{args.out_prefix}.gene.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
