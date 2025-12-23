BASE=/data/holwm/home/gaolong/projects/RiboSuite/reference/human/GRCh38
GTF=${BASE}/gencode.v45.annotation.gtf
FASTA=${BASE}/GRCh38.primary_assembly.genome.fa

OUTDIR=${BASE}/contam
mkdir -p ${OUTDIR}

# # Biotypes to include (edit here if you want)
# BIOTYPES='rRNA|Mt_rRNA|tRNA|Mt_tRNA|snRNA|snoRNA|scaRNA|vault_RNA|Y_RNA|misc_RNA|rRNA_pseudogene|tRNA_pseudogene'

# # Make a BED of exons from these biotypes (0-based BED coordinates)
# awk -v pat="^("${BIOTYPES}")$" '
# BEGIN{FS=OFS="\t"}
# $3=="exon" {
#   gene_type=""; gene_id=""; gene_name=""; transcript_id="";
#   if (match($0,/gene_type "([^"]+)"/,a)) gene_type=a[1];
#   if (gene_type ~ pat) {
#     if (match($0,/gene_id "([^"]+)"/,b)) gene_id=b[1];
#     if (match($0,/gene_name "([^"]+)"/,c)) gene_name=c[1];
#     if (match($0,/transcript_id "([^"]+)"/,d)) transcript_id=d[1];
#     # BED: chrom, start(0-based), end, name, score, strand
#     name=gene_type"|"gene_name"|"gene_id"|"transcript_id"|exon";
#     print $1, $4-1, $5, name, 0, $7;
#   }
# }
# ' ${GTF} \
# | sort -k1,1 -k2,2n \
# > ${OUTDIR}/human_contam.exons.bed



# bedtools getfasta \
#   -fi ${FASTA} \
#   -bed ${OUTDIR}/human_contam.exons.bed \
#   -s \
#   -name \
#   -fo ${OUTDIR}/human_contam.gencode_v45.fa


IDXDIR=${OUTDIR}/bowtie2_idx
mkdir -p ${IDXDIR}

bowtie2-build \
  ${OUTDIR}/human_contam.gencode_v45.fa \
  ${IDXDIR}/human_contam