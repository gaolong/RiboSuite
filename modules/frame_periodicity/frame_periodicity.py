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
        description="Fast start-window frame periodicity (protein-coding, longest CDS per gene)"
    )
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--sample", required=True)

    ap.add_argument("--start_window", default="-30,10",
                    help="Start window relative to CDS start (nt), e.g. -30,10")
    ap.add_argument("--min_cds_len", type=int, default=150,
                    help="Minimum CDS length to keep transcript")
    ap.add_argument("--min_start_reads", type=int, default=0,
                    help="Minimum reads in start window to keep gene (0 disables)")

    return ap.parse_args()

# --------------------------------------------------
# Offsets
# --------------------------------------------------
def load_offsets(path):
    df = pd.read_csv(path, sep="\t")
    col = "psite_offset" if "psite_offset" in df.columns else "offset"
    return dict(zip(df["read_length"], df[col]))

# --------------------------------------------------
# GTF parsing
# --------------------------------------------------
def parse_attrs(attr):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr))

def load_longest_cds_per_gene(gtf, min_cds_len):
    """
    Returns:
      cds_models: list of dicts:
        {
          chrom, strand, gene_id,
          cds_exons: [(start,end,phase)],
          cds_len,
          cds_start
        }
    """
    tx_cds = defaultdict(list)
    tx_gene = {}
    tx_gene_type = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, phase, attr = line.rstrip().split("\t")
            attrs = parse_attrs(attr)

            if feature == "gene":
                tx_gene_type[attrs["gene_id"]] = attrs.get("gene_type") or attrs.get("gene_biotype")

            if feature != "CDS":
                continue

            gene = attrs.get("gene_id")
            tx = attrs.get("transcript_id")
            if gene is None or tx is None:
                continue

            start0 = int(start) - 1
            end0 = int(end)
            ph = 0 if phase == "." else int(phase)

            tx_cds[(gene, tx, chrom, strand)].append((start0, end0, ph))
            tx_gene[tx] = gene

    # choose longest CDS transcript per gene
    best = {}

    for (gene, tx, chrom, strand), exons in tx_cds.items():
        if tx_gene_type.get(gene) != "protein_coding":
            continue

        exons.sort(key=lambda x: x[0])
        cds_len = sum(e - s for s, e, _ in exons)
        if cds_len < min_cds_len:
            continue

        if gene not in best or cds_len > best[gene]["cds_len"]:
            cds_start = exons[0][0] if strand == "+" else exons[-1][1] - 1
            best[gene] = {
                "chrom": chrom,
                "strand": strand,
                "gene_id": gene,
                "cds_exons": exons,
                "cds_len": cds_len,
                "cds_start": cds_start
            }

    return list(best.values())

# --------------------------------------------------
# Build spliced start-window segments
# --------------------------------------------------
def build_start_window_segments(cds_exons, cds_start, strand, wL, wR):
    """
    Map transcript-relative start window onto genomic CDS exons.
    Returns list of (start,end,phase).
    """
    segs = []

    if strand == "+":
        cursor = cds_start + wL
        remaining = wR - wL + 1
        for s, e, ph in cds_exons:
            if cursor >= e:
                continue
            if cursor < s:
                cursor = s
            take = min(e - cursor, remaining)
            if take > 0:
                segs.append((cursor, cursor + take, ph))
                remaining -= take
                cursor += take
            if remaining <= 0:
                break
    else:
        cursor = cds_start - wL
        remaining = wR - wL + 1
        for s, e, ph in reversed(cds_exons):
            if cursor < s:
                continue
            if cursor >= e:
                cursor = e - 1
            take = min(cursor - s + 1, remaining)
            if take > 0:
                segs.append((cursor - take + 1, cursor + 1, ph))
                remaining -= take
                cursor -= take
            if remaining <= 0:
                break

    return segs

# --------------------------------------------------
# Frame calculation
# --------------------------------------------------
def frame_with_phase(psite, s, e, ph, strand):
    if strand == "+":
        return (psite - s + ph) % 3
    else:
        return ((e - 1 - psite) + ph) % 3

# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    args = parse_args()
    offsets = load_offsets(args.offset)
    wL, wR = map(int, args.start_window.split(","))

    cds_models = load_longest_cds_per_gene(args.gtf, args.min_cds_len)

    bam = pysam.AlignmentFile(args.bam, "rb")

    counts = defaultdict(int)

    for m in cds_models:
        segs = build_start_window_segments(
            m["cds_exons"], m["cds_start"], m["strand"], wL, wR
        )

        gene_reads = 0

        for s, e, ph in segs:
            for read in bam.fetch(m["chrom"], s, e):
                if read.is_unmapped:
                    continue
                rl = read.query_length
                if rl not in offsets:
                    continue

                offset = offsets[rl]
                psite = (
                    read.reference_start + offset
                    if m["strand"] == "+"
                    else read.reference_end - offset - 1
                )

                if not (s <= psite < e):
                    continue

                frame = frame_with_phase(psite, s, e, ph, m["strand"])
                counts[(rl, frame)] += 1
                gene_reads += 1

        if args.min_start_reads and gene_reads < args.min_start_reads:
            continue

    # output
    rows = [
        {"sample": args.sample, "read_length": rl, "frame": fr, "count": c}
        for (rl, fr), c in counts.items()
    ]

    df = pd.DataFrame(rows).sort_values(["read_length", "frame"])
    df["fraction"] = df["count"] / df.groupby("read_length")["count"].transform("sum")

    df.to_csv(f"{args.sample}.periodicity.by_region_by_length.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
