import pysam
import argparse
import collections
import sys

def parse_args():
    p = argparse.ArgumentParser()
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


def load_start_codons(gtf):
    """
    Parse GTF and return a list of collapsed translation start sites.

    Each returned entry corresponds to ONE unique genomic start site
    (chrom, start_pos, strand), with all supporting transcripts combined.

    Returns:
        List of tuples:
        (chrom, start_pos, strand, gene_names, transcript_ids)

        where gene_names and transcript_ids are comma-separated strings.
    """

    # Step 1: collect CDSs per transcript
    cds_by_tx = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue

            if fields[2] != "CDS":
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # GTF is 1-based
            end = int(fields[4]) - 1
            strand = fields[6]

            # parse attributes into dict
            attrs = {}
            for item in fields[8].split(";"):
                item = item.strip()
                if not item:
                    continue
                key, val = item.split(" ", 1)
                attrs[key] = val.strip('"')

            tx_id = attrs.get("transcript_id")
            gene_name = attrs.get("gene_name", "NA")

            if tx_id is None:
                continue

            if tx_id not in cds_by_tx:
                cds_by_tx[tx_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_name": gene_name,
                    "cds": []
                }

            cds_by_tx[tx_id]["cds"].append((start, end))

    # Step 2: determine first CDS per transcript
    # and collapse by (chrom, start_pos, strand)
    start_dict = {}  # key = (chrom, start_pos, strand)

    for tx_id, info in cds_by_tx.items():
        chrom = info["chrom"]
        strand = info["strand"]
        gene_name = info["gene_name"]
        cds_list = info["cds"]

        if strand == "+":
            start_pos = min(cds_list, key=lambda x: x[0])[0]
        else:
            start_pos = max(cds_list, key=lambda x: x[1])[1]

        key = (chrom, start_pos, strand)

        if key not in start_dict:
            start_dict[key] = {
                "chrom": chrom,
                "start_pos": start_pos,
                "strand": strand,
                "gene_names": set(),
                "transcript_ids": set()
            }

        start_dict[key]["gene_names"].add(gene_name)
        start_dict[key]["transcript_ids"].add(tx_id)

    # Step 3: format output
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

    # distance counts per read length
    dist = collections.defaultdict(list)

    for chrom, start_pos, strand in starts:
        for read in bam.fetch(
            chrom,
            start_pos - args.window_up,
            start_pos + args.window_down
        ):
            if read.is_unmapped:
                continue

            read_len = read.query_length
            if read_len < args.min_len or read_len > args.max_len:
                continue

            if strand == "+":
                read_5p = read.reference_start
                rel = read_5p - start_pos
            else:
                read_5p = read.reference_end - 1
                rel = start_pos - read_5p

            if -args.window_up <= rel <= args.window_down:
                dist[read_len].append(-rel)

    with open(args.out, "w") as out:
        out.write("sample_id\tread_length\tpsite_offset\tn_reads\n")
        for L in sorted(dist):
            if len(dist[L]) < args.min_reads:
                continue
            counter = collections.Counter(dist[L])
            offset, n = counter.most_common(1)[0]
            out.write(
                f"{args.sample_id}\t{L}\t{offset}\t{len(dist[L])}\n"
            )

    bam.close()


if __name__ == "__main__":
    main()