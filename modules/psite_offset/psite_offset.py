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
    Return list of (chrom, pos, strand) for CDS start codons
    """
    starts = []
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "CDS":
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4]) - 1
            strand = fields[6]
            if strand == "+":
                starts.append((chrom, start, strand))
            else:
                starts.append((chrom, end, strand))
    return starts

def main():
    args = parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    starts = load_start_codons(args.gtf)

    # distance counts per read length
    dist = collections.defaultdict(list)

    for chrom, start_pos, strand in starts:
        for read in bam.fetch(chrom, start_pos - args.window_up, start_pos + args.window_down):
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
            out.write(f"{args.sample_id}\t{L}\t{offset}\t{len(dist[L])}\n")

    bam.close()

if __name__ == "__main__":
    main()