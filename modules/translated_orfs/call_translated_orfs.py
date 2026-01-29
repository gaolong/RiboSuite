#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
from collections import defaultdict
import re

# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Rule-based identification of translated ORFs "
                    "using start-codon–anchored CDS frame periodicity"
    )
    p.add_argument("--bam", required=True, help="Deduplicated, sorted BAM")
    p.add_argument("--gtf", required=True, help="GTF annotation")
    p.add_argument("--psite_offsets", required=True, help="P-site offset table")
    p.add_argument("--out_prefix", required=True)

    p.add_argument("--min_psites", type=int, default=30)
    p.add_argument("--min_len_codons", type=int, default=30)
    p.add_argument("--min_frame_fraction", type=float, default=0.7)
    p.add_argument("--mapq", type=int, default=100)

    return p.parse_args()

# --------------------------------------------------
# offsets
# --------------------------------------------------
def load_psite_offsets(path):
    df = pd.read_csv(path, sep="\t")
    col = "psite_offset" if "psite_offset" in df.columns else "offset"
    return dict(zip(df["read_length"].astype(int), df[col].astype(int)))

# --------------------------------------------------
# GTF helpers
# --------------------------------------------------
def parse_attrs(attr):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr))

# --------------------------------------------------
# genomic → spliced transcript coordinate
# --------------------------------------------------
def genomic_to_tx_spliced(pos, exons, strand):
    for s, e, tx_s in exons:
        if s <= pos < e:
            return tx_s + (pos - s)
    return None

# --------------------------------------------------
# load transcript models (adapted from frame_periodicity.py)
# --------------------------------------------------
def load_tx_models(gtf, min_len_codons):
    tx_exons = defaultdict(list)
    tx_cds = defaultdict(list)
    tx_start = defaultdict(list)
    tx_stop = defaultdict(list)
    gene_type = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, phase, attr = (
                line.rstrip().split("\t")
            )
            attrs = parse_attrs(attr)

            if feature == "gene":
                gene_type[attrs["gene_id"]] = (
                    attrs.get("gene_type") or attrs.get("gene_biotype")
                )
                continue

            gene = attrs.get("gene_id")
            tx = attrs.get("transcript_id")
            if not gene or not tx:
                continue
            if gene_type.get(gene) != "protein_coding":
                continue

            start0 = int(start) - 1
            end0 = int(end)
            key = (gene, tx, chrom, strand)

            if feature == "exon":
                tx_exons[key].append((start0, end0))
            elif feature == "CDS":
                ph = 0 if phase == "." else int(phase)
                tx_cds[key].append((start0, end0, ph))
            elif feature == "start_codon":
                tx_start[key].append((start0, end0))
            elif feature == "stop_codon":
                tx_stop[key].append((start0, end0))

    models = {}

    for key, exons in tx_exons.items():
        gene, tx, chrom, strand = key

        if key not in tx_cds or key not in tx_start or key not in tx_stop:
            continue

        cds_len = sum(e - s - ph for s, e, ph in tx_cds[key])
        if cds_len < min_len_codons * 3:
            continue

        # sort exons
        exons.sort(key=lambda x: x[0])

        # build spliced transcript coordinates
        tx_pos = []
        cursor = 0
        for s, e in exons:
            tx_pos.append((s, e, cursor))
            cursor += (e - s)

        # start codon
        if strand == "+":
            sc_g = min(s for s, e in tx_start[key])
        else:
            sc_g = max(e for s, e in tx_start[key]) - 1

        sc_tx = genomic_to_tx_spliced(sc_g, tx_pos, strand)
        if sc_tx is None:
            continue

        # stop codon
        if strand == "+":
            st_g = max(e for s, e in tx_stop[key]) - 1
        else:
            st_g = min(s for s, e in tx_stop[key])

        st_tx = genomic_to_tx_spliced(st_g, tx_pos, strand)
        if st_tx is None:
            continue

        cds_lo_tx = min(sc_tx, st_tx)
        cds_hi_tx = max(sc_tx, st_tx)

        models[tx] = {
            "gene": gene,
            "tx": tx,
            "chrom": chrom,
            "strand": strand,
            "exons": tx_pos,
            "start_codon_tx": sc_tx,
            "cds_lo_tx": cds_lo_tx,
            "cds_hi_tx": cds_hi_tx,
            "fetch_start": min(s for s, e in exons),
            "fetch_end": max(e for s, e in exons),
        }

    return models

# --------------------------------------------------
# core logic
# --------------------------------------------------
def assign_psites_to_tx(bam, offsets, m, mapq):
    fcounts = [0, 0, 0]
    total = 0
    strand = m["strand"]
    sc_tx = m["start_codon_tx"]

    max_offset = max(offsets.values())
    fs = max(0, m["fetch_start"] - max_offset)
    fe = m["fetch_end"] + max_offset

    for read in bam.fetch(m["chrom"], fs, fe):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < mapq:
            continue

        rl = read.query_length
        if rl not in offsets:
            continue

        psite = (
            read.reference_start + offsets[rl]
            if strand == "+"
            else read.reference_end - offsets[rl] - 1
        )

        tx_pos = genomic_to_tx_spliced(psite, m["exons"], strand)
        if tx_pos is None:
            continue

        if not (m["cds_lo_tx"] <= tx_pos < m["cds_hi_tx"]):
            continue

        # FRAME (identical to frame_periodicity.py)
        if strand == "+":
            frame = (tx_pos - sc_tx) % 3
        else:
            frame = (sc_tx - tx_pos) % 3

        fcounts[frame] += 1
        total += 1

    return total, fcounts

# --------------------------------------------------
# rules
# --------------------------------------------------
def rule_min_psites(total, min_psites):
    return total >= min_psites


def rule_frame_dominance(fcounts, min_fraction):
    s = sum(fcounts)
    if s == 0:
        return False, 0.0
    top = max(fcounts)
    frac = top / s
    return frac >= min_fraction, frac

# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()

    offsets = load_psite_offsets(args.psite_offsets)
    bam = pysam.AlignmentFile(args.bam, "rb")

    tx_models = load_tx_models(args.gtf, args.min_len_codons)

    rows = []

    for tx, m in tx_models.items():
        total, fcounts = assign_psites_to_tx(
            bam, offsets, m, args.mapq
        )

        rules_passed = []
        rules_failed = []

        if rule_min_psites(total, args.min_psites):
            rules_passed.append("min_psites")
        else:
            rules_failed.append("min_psites")

        ok_frame, frac = rule_frame_dominance(
            fcounts, args.min_frame_fraction
        )
        if ok_frame:
            rules_passed.append("frame_dominance")
        else:
            rules_failed.append("frame_dominance")

        rows.append({
            "gene_id": m["gene"],
            "transcript_id": tx,
            "psite_total": total,
            "psite_f0": fcounts[0],
            "psite_f1": fcounts[1],
            "psite_f2": fcounts[2],
            "frame_fraction": round(frac, 3),
            "rules_passed": ",".join(rules_passed),
            "rules_failed": ",".join(rules_failed),
        })

    df = pd.DataFrame(rows)
    df.to_csv(
        f"{args.out_prefix}.translated_orfs.tsv",
        sep="\t",
        index=False
    )

    print(f"[translated_orfs] written {args.out_prefix}.translated_orfs.tsv")

if __name__ == "__main__":
    main()
