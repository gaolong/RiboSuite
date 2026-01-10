#!/usr/bin/env python3

import pysam
import argparse
import pandas as pd
from collections import defaultdict
import re

# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Frame periodicity QC (riboWaltz / plastid logic)"
    )
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min_cds_len", type=int, default=150)
    ap.add_argument("--top_tx", type=int, default=1000)
    ap.add_argument("--mapq", type=int, default=10)
    return ap.parse_args()

# --------------------------------------------------
# offsets
# --------------------------------------------------
def load_offsets(path):
    df = pd.read_csv(path, sep="\t")
    col = "psite_offset" if "psite_offset" in df.columns else "offset"
    return dict(zip(df["read_length"].astype(int), df[col].astype(int)))

# --------------------------------------------------
# GTF parsing
# --------------------------------------------------
def parse_attrs(attr):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr))

def load_tx_models(gtf, min_cds_len):
    tx_cds = defaultdict(list)
    gene_type = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, phase, attr = line.rstrip().split("\t")
            attrs = parse_attrs(attr)

            if feature == "gene":
                gene_type[attrs["gene_id"]] = attrs.get("gene_type") or attrs.get("gene_biotype")

            if feature != "CDS":
                continue

            gene = attrs.get("gene_id")
            tx   = attrs.get("transcript_id")
            if not gene or not tx:
                continue

            if gene_type.get(gene) != "protein_coding":
                continue

            start0 = int(start) - 1
            end0   = int(end)
            ph = 0 if phase == "." else int(phase)

            tx_cds[(gene, tx, chrom, strand)].append((start0, end0, ph))

    models = []

    for (gene, tx, chrom, strand), exons in tx_cds.items():
        exons.sort(key=lambda x: x[0])

        cds_len = sum(e - s for s, e, _ in exons)
        if cds_len < min_cds_len:
            continue

        # build transcript coordinate map
        tx_pos = []
        cursor = 0
        for s, e, ph in exons:
            length = e - s
            tx_pos.append((s, e, cursor, ph))
            cursor += length

        models.append({
            "gene": gene,
            "tx": tx,
            "chrom": chrom,
            "strand": strand,
            "exons": tx_pos,
            "cds_len": cds_len
        })

    return models

# --------------------------------------------------
# genomic → transcript coordinate
# --------------------------------------------------
def genomic_to_tx(psite, exons, strand):
    """
    Convert genomic P-site to transcript coordinate
    """
    if strand == "+":
        for s, e, tx_start, _ in exons:
            if s <= psite < e:
                return tx_start + (psite - s)
    else:
        for s, e, tx_start, _ in reversed(exons):
            if s <= psite < e:
                return tx_start + (e - psite - 1)
    return None

# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()
    offsets = load_offsets(args.offset)
    bam = pysam.AlignmentFile(args.bam, "rb")

    models = load_tx_models(args.gtf, args.min_cds_len)

    # --------------------------------------------------
    # PASS 1: select top transcripts
    # --------------------------------------------------
    tx_counts = defaultdict(int)

    for m in models:
        chrom, strand = m["chrom"], m["strand"]
        tx_id = (m["gene"], m["tx"])

        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq:
                continue

            rl = read.query_length
            if rl not in offsets:
                continue

            psite = (
                read.reference_start + offsets[rl]
                if strand == "+"
                else read.reference_end - offsets[rl] - 1
            )

            tx_pos = genomic_to_tx(psite, m["exons"], strand)
            if tx_pos is not None:
                tx_counts[tx_id] += 1

    top_tx = set(
        tx for tx, _ in sorted(
            tx_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:args.top_tx]
    )

    # --------------------------------------------------
    # PASS 2: frame periodicity
    # --------------------------------------------------
    counts = defaultdict(int)

    for m in models:
        tx_id = (m["gene"], m["tx"])
        if tx_id not in top_tx:
            continue

        chrom, strand = m["chrom"], m["strand"]

        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq:
                continue

            rl = read.query_length
            if rl not in offsets:
                continue

            psite = (
                read.reference_start + offsets[rl]
                if strand == "+"
                else read.reference_end - offsets[rl] - 1
            )

            tx_pos = genomic_to_tx(psite, m["exons"], strand)
            if tx_pos is None:
                continue

            frame = tx_pos % 3
            counts[(rl, frame)] += 1

    # --------------------------------------------------
    # output
    # --------------------------------------------------
    rows = [
        dict(
            sample=args.sample,
            read_length=rl,
            frame=fr,
            count=c
        )
        for (rl, fr), c in counts.items()
    ]

    df = pd.DataFrame(rows)
    df["fraction"] = df["count"] / df.groupby("read_length")["count"].transform("sum")

    out = f"{args.sample}.periodicity.ribowaltz_style.tsv"
    df.to_csv(out, sep="\t", index=False)

    print(f"[periodicity] written {out}")

if __name__ == "__main__":
    main()
