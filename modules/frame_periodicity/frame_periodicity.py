#!/usr/bin/env python3

import pysam
import argparse
import pandas as pd
from collections import defaultdict
import re

# --- plotting (HPC / headless safe) ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, FixedFormatter


# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Frame periodicity QC for 5'UTR, CDS and 3'UTR "
                    "(CDS-anchored frames, optional phase support)"
    )
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min_cds_len", type=int, default=0)
    ap.add_argument("--top_tx", type=int, default=15000)
    ap.add_argument("--mapq", type=int, default=0)

    # NEW: optional phase usage
    ap.add_argument(
        "--use_phase",
        action="store_true",
        help="Use CDS phase from GTF for frame calculation (ribowaltz-style)"
    )

    return ap.parse_args()


# --------------------------------------------------
# offsets
# --------------------------------------------------
def load_offsets(path):
    df = pd.read_csv(path, sep="\t")
    col = "psite_offset" if "psite_offset" in df.columns else "offset"
    return dict(zip(df["read_length"].astype(int), df[col].astype(int)))


# --------------------------------------------------
# intervals
# --------------------------------------------------
def in_interval(pos, intervals):
    for s, e in intervals:
        if s <= pos < e:
            return True
    return False


# --------------------------------------------------
# GTF helpers
# --------------------------------------------------
def parse_attrs(attr):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr))


def load_start_codons(gtf):
    start_tx = set()
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            _, _, feature, *_ , attr = line.rstrip().split("\t")
            if feature != "start_codon":
                continue
            attrs = parse_attrs(attr)
            if "gene_id" in attrs and "transcript_id" in attrs:
                start_tx.add((attrs["gene_id"], attrs["transcript_id"]))
    return start_tx


# --------------------------------------------------
# genomic → spliced transcript coordinate
# --------------------------------------------------
def genomic_to_tx_spliced(pos, exons, strand):
    for s, e, tx_s in exons:
        if s <= pos < e:
            return tx_s + (pos - s)
    return None


# --------------------------------------------------
# load transcript models
# --------------------------------------------------
def load_tx_models(gtf, min_cds_len):
    """
    Load transcript models for ribosome profiling.

    Keeps ONLY transcripts that have:
      - at least one CDS
      - start_codon
      - stop_codon

    CRITICAL DEFINITIONS
    --------------------
    - Frame 0 is anchored at START CODON (not CDS start)
    - Ribo-seq CDS region is [start_codon, stop_codon)
    - Annotated CDS bounds are kept for QC ONLY

    Returns
    -------
    models  : list of transcript models for read counting
    tx_meta : per-transcript metadata
    """

    tx_exons = defaultdict(list)
    tx_cds = defaultdict(list)
    tx_utrs = defaultdict(list)
    tx_start_codon = defaultdict(list)
    tx_stop_codon = defaultdict(list)
    gene_type = {}

    # --------------------------------------------------
    # Parse GTF
    # --------------------------------------------------
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

            elif feature == "UTR":
                tx_utrs[key].append((start0, end0))

            elif feature == "start_codon":
                tx_start_codon[key].append((start0, end0))

            elif feature == "stop_codon":
                tx_stop_codon[key].append((start0, end0))
    
    models = []
    tx_meta = {}

    # --------------------------------------------------
    # Build transcript models
    # --------------------------------------------------
    for key, exons in tx_exons.items():
        gene, tx, chrom, strand = key

        if key not in tx_cds:
            continue
        if key not in tx_start_codon:
            continue
        if key not in tx_stop_codon:
            continue

        # --------------------------------------------------
        # Sort exons in genomic order
        # --------------------------------------------------
        exons.sort(key=lambda x: x[0])

        # --------------------------------------------------
        # Annotated CDS length (QC only)
        # --------------------------------------------------
        cds_intervals = tx_cds[key]
        cds_len = sum(e - s - ph for s, e, ph in cds_intervals)
        if cds_len < min_cds_len:
            continue

        # --------------------------------------------------
        # Build spliced transcript coordinates
        # --------------------------------------------------
        tx_pos = []
        cursor = 0
        for s, e in exons:
            tx_pos.append((s, e, cursor))
            cursor += (e - s)

        # --------------------------------------------------
        # Resolve START codon (genomic → tx)
        # --------------------------------------------------
        sc_intervals = tx_start_codon[key]
        if strand == "+":
            start_codon_g = min(s for s, e in sc_intervals)
        else:
            start_codon_g = max(e for s, e in sc_intervals) - 1

        start_codon_tx = genomic_to_tx_spliced(
            start_codon_g, tx_pos, strand
        )
        if start_codon_tx is None:
            continue

        # --------------------------------------------------
        # Resolve STOP codon (genomic → tx)
        # --------------------------------------------------
        stop_intervals = tx_stop_codon[key]
        if strand == "+":
            stop_codon_g = max(e for s, e in stop_intervals) - 1
        else:
            stop_codon_g = min(s for s, e in stop_intervals)

        stop_codon_tx = genomic_to_tx_spliced(
            stop_codon_g, tx_pos, strand
        )
        if stop_codon_tx is None:
            continue

        # --------------------------------------------------
        # Translation bounds (USED FOR RIBO-SEQ)
        # --------------------------------------------------
        cds_lo_tx = min(start_codon_tx, stop_codon_tx)
        cds_hi_tx = max(start_codon_tx, stop_codon_tx)

        # --------------------------------------------------
        # Annotated CDS bounds (QC ONLY)
        # --------------------------------------------------
        cds_start_g = min(s + ph for s, e, ph in cds_intervals)
        cds_end_g = max(e for s, e, ph in cds_intervals)

        annot_cds_start_tx = genomic_to_tx_spliced(
            cds_start_g, tx_pos, strand
        )
        annot_cds_end_tx_pos = genomic_to_tx_spliced(
            cds_end_g - 1, tx_pos, strand
        )

        annot_cds_end_tx = (
            annot_cds_end_tx_pos + 1
            if annot_cds_end_tx_pos is not None
            else None
        )

        # --------------------------------------------------
        # UTRs → transcript coordinates
        # --------------------------------------------------
        utr_tx = []
        for s, e in tx_utrs.get(key, []):
            tx_s = genomic_to_tx_spliced(s, tx_pos, strand)
            tx_e = genomic_to_tx_spliced(e - 1, tx_pos, strand)
            if tx_s is not None and tx_e is not None:
                utr_tx.append((tx_s, tx_e + 1))

        # --------------------------------------------------
        # Store model (COUNTING)
        # --------------------------------------------------
        models.append({
            "gene": gene,
            "tx": tx,
            "chrom": chrom,
            "strand": strand,
            "exons": tx_pos,

            # Translation anchor (CRITICAL)
            "start_codon_tx": start_codon_tx,
            "stop_codon_tx": stop_codon_tx,

            # Ribo-seq CDS bounds (USED)
            "cds_lo_tx": cds_lo_tx,
            "cds_hi_tx": cds_hi_tx,

            # Annotation-only CDS bounds (DO NOT USE FOR FRAME)
            "annot_cds_start_tx": annot_cds_start_tx,
            "annot_cds_end_tx": annot_cds_end_tx,

            "utr_tx": utr_tx,
            "fetch_start": min(s for s, e in exons),
            "fetch_end": max(e for s, e in exons),
        })

        # --------------------------------------------------
        # Metadata
        # --------------------------------------------------
        tx_meta[tx] = {
            "gene": gene,
            "chrom": chrom,
            "strand": strand,

            "start_codon_genomic": start_codon_g,
            "start_codon_tx": start_codon_tx,
            "stop_codon_genomic": stop_codon_g,
            "stop_codon_tx": stop_codon_tx,

            "cds_lo_tx": cds_lo_tx,
            "cds_hi_tx": cds_hi_tx,

            "annot_cds_start_tx": annot_cds_start_tx,
            "annot_cds_end_tx": annot_cds_end_tx,

            "cds_len": cds_len,
            "tx_len": cursor,
            "fetch_start": min(s for s, e in exons),
            "fetch_end": max(e for s, e in exons),
        }

    return models


# --------------------------------------------------
# main
# --------------------------------------------------
def main():
    args = parse_args()

    offsets = load_offsets(args.offset)
    max_offset = max(offsets.values())
    bam = pysam.AlignmentFile(args.bam, "rb")

    models = load_tx_models(args.gtf, args.min_cds_len)
    start_tx = load_start_codons(args.gtf)
    models = [m for m in models if (m["gene"], m["tx"]) in start_tx]

    if not models:
        raise RuntimeError("No valid transcripts found")

    # --------------------------------------------------
    # PASS 1: select top transcripts by CDS coverage
    # --------------------------------------------------
    tx_counts = defaultdict(int)

    for m in models:
        chrom, strand = m["chrom"], m["strand"]
        tx_id = (m["gene"], m["tx"])

        fs = max(0, m["fetch_start"] - max_offset)
        fe = m["fetch_end"] + max_offset

        for read in bam.fetch(chrom, fs, fe):
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

            tx_pos = genomic_to_tx_spliced(psite, m["exons"], strand)
            if tx_pos is None:
                continue

            if m["cds_lo_tx"] <= tx_pos < m["cds_hi_tx"]:
                tx_counts[tx_id] += 1

    top_tx = {
        tx for tx, _ in sorted(
            tx_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:args.top_tx]
    }

    # --------------------------------------------------
    # PASS 2: frame periodicity
    # --------------------------------------------------
    counts = defaultdict(int)

    for m in models:
        if (m["gene"], m["tx"]) not in top_tx:
            continue

        chrom, strand = m["chrom"], m["strand"]
        fs = max(0, m["fetch_start"] - max_offset)
        fe = m["fetch_end"] + max_offset
        start_codon_tx = m["start_codon_tx"]

        for read in bam.fetch(chrom, fs, fe):
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

            tx_pos = genomic_to_tx_spliced(psite, m["exons"], strand)
            if tx_pos is None:
                continue

            in_cds = m["cds_lo_tx"] <= tx_pos < m["cds_hi_tx"]
            in_utr = in_interval(tx_pos, m["utr_tx"])

            if not (in_cds or in_utr):
                continue

            # FRAME CALCULATION
            if strand == "+":
                frame = (tx_pos - start_codon_tx) % 3
            else:
                frame = (start_codon_tx - tx_pos) % 3

            if in_cds:
                region = "CDS"
            elif tx_pos < m["cds_lo_tx"]:
                region = "5UTR"
            else:
                region = "3UTR"

            counts[(region, rl, frame)] += 1

    # --------------------------------------------------
    # output
    # --------------------------------------------------
    df = pd.DataFrame([
        dict(sample=args.sample, region=r, read_length=rl, frame=f, count=c)
        for (r, rl, f), c in counts.items()
    ])
    
    out_tsv = f"{args.sample}.periodicity.by_region_by_length.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[periodicity] written {out_tsv}")

    if df.empty:
        return

    df["fraction"] = df["count"] / df.groupby(
        ["region", "read_length"]
    )["count"].transform("sum")

    # --------------------------------------------------
    # plotting
    # --------------------------------------------------
    regions = ["5UTR", "CDS", "3UTR"]
    titles = ["5′ UTR", "CDS", "3′ UTR"]

    cmap = LinearSegmentedColormap.from_list("white_red", ["#ffffff", "#ff0000"])
    lengths = sorted(df["read_length"].unique())

    fig, axes = plt.subplots(
        ncols=3,
        figsize=(11, max(4, len(lengths) * 0.35)),
        constrained_layout=True
    )

    for ax, region, title in zip(axes, regions, titles):
        sub = df[df["region"] == region]

        mat = (
            sub.pivot_table(
                index="read_length",
                columns="frame",
                values="fraction",
                fill_value=0.0
            )
            .reindex(lengths, fill_value=0.0)
        )

        im = ax.imshow(
            mat.values,
            aspect="auto",
            cmap=cmap,
            vmin=0,
            vmax=1,
            origin="upper"
        )

        ax.set_title(title)

        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["Frame 0", "Frame 1", "Frame 2"])

        ax.set_yticks(range(len(lengths)))
        if ax is axes[0]:
            ax.set_yticklabels(lengths)
            ax.set_ylabel("Read length")
        else:
            ax.set_yticklabels([])

    fig.colorbar(im, ax=axes, fraction=0.035, pad=0.02).set_label("Fraction")
    fig.suptitle(f"{args.sample} frame periodicity")

    out_png = f"{args.sample}.periodicity.by_region_by_length.heatmap.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"[periodicity] written {out_png}")


if __name__ == "__main__":
    main()
