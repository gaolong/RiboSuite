#!/usr/bin/env python

import pysam
import argparse
import collections


def parse_args():
    p = argparse.ArgumentParser(
        description="Estimate P-site offsets using CDS start sites (mode-based)"
    )
    p.add_argument("--bam", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--sample_id", required=True)
    p.add_argument("--window_up", type=int, default=50)
    p.add_argument("--window_down", type=int, default=20)
    p.add_argument("--min_len", type=int, default=26)
    p.add_argument("--max_len", type=int, default=34)
    p.add_argument("--min_reads", type=int, default=100)
    p.add_argument("--out", required=True)
    return p.parse_args()


def parse_gtf_attributes(attr_str):
    attrs = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if not item:
            continue
        key, val = item.split(" ", 1)
        attrs[key] = val.strip('"')
    return attrs


def load_start_codons(gtf):
    """
    Updated behavior:
    - use explicit start_codon features
    - collapse transcripts sharing the same genomic start
    """

    start_dict = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            if feature != "start_codon":
                continue

            chrom = fields[0]
            strand = fields[6]

            # GTF is 1-based, inclusive
            start = int(fields[3]) - 1
            end = int(fields[4]) - 1

            attrs = parse_gtf_attributes(fields[8])
            tx_id = attrs.get("transcript_id")
            gene_name = attrs.get("gene_name", "NA")

            if tx_id is None:
                continue

            # choose the correct genomic coordinate for CDS start
            if strand == "+":
                start_pos = start
            else:
                start_pos = end

            key = (chrom, start_pos, strand)

            start_dict.setdefault(key, {
                "chrom": chrom,
                "start_pos": start_pos,
                "strand": strand,
                "gene_names": set(),
                "transcript_ids": set()
            })

            start_dict[key]["gene_names"].add(gene_name)
            start_dict[key]["transcript_ids"].add(tx_id)

    starts = []
    for v in start_dict.values():
        starts.append((
            v["chrom"],
            v["start_pos"],
            v["strand"],
            ",".join(sorted(v["gene_names"])),
            ",".join(sorted(v["transcript_ids"]))
        ))

    return starts



def main():
    args = parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    starts = load_start_codons(args.gtf)

    # store distance = (CDS_start - read_5p)
    dist = collections.defaultdict(list)

    for chrom, start_pos, strand, _, _ in starts:
        for read in bam.fetch(
            chrom,
            start_pos - args.window_up,
            start_pos + args.window_down
        ):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            L = read.query_length
            if L < args.min_len or L > args.max_len:
                continue

            # correct 5′ end of the read
            read_5p = read.reference_end - 1 if read.is_reverse else read.reference_start

            # enforce strand consistency
            if read.is_reverse != (strand == "-"):
                continue

            # distance from CDS start to read 5′ end (transcript direction)
            if strand == "+":
                rel = read_5p - start_pos
            else:
                rel = start_pos - read_5p

            if -args.window_up <= rel <= args.window_down:
                dist[L].append(-rel)  # IMPORTANT: keep sign convention

    with open(args.out, "w") as out:
        out.write("sample_id\tread_length\tpsite_offset\tn_reads\n")

        for L in sorted(dist):
            if len(dist[L]) < args.min_reads:
                continue

            counter = collections.Counter(dist[L])
            offset = counter.most_common(1)[0][0]

            # compute frame-0 fraction using chosen offset
            f0 = f1 = f2 = 0
            for d in dist[L]:
                # P-site relative to CDS start = offset - d
                frame = (offset - d) % 3
                if frame == 0:
                    f0 += 1
                elif frame == 1:
                    f1 += 1
                else:
                    f2 += 1

            total = f0 + f1 + f2
            frame0_frac = f0 / total if total > 0 else 0.0

            out.write(
                f"{args.sample_id}\t{L}\t{offset}\t{len(dist[L])}\n"
            )

    bam.close()


if __name__ == "__main__":
    main()
