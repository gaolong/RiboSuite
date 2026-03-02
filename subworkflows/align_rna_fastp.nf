nextflow.enable.dsl = 2

include { FASTP }          from '../modules/fastp/main.nf'
include { STAR_RNA_ALIGN } from '../modules/star_rna/main.nf'

/*
 * RNA-seq preprocess + alignment
 *
 * Input:
 *   reads_ch: tuple(meta, read1, read2)   (paired-end)
 *   star_index: path to STAR genome index
 *
 * Output:
 *   bam:        tuple(meta, bam)
 *   sj:         path SJ.out.tab
 *   gene_counts: path ReadsPerGene.out.tab
 *   log:        path Log.final.out
 *   fastp_json/html
 */
workflow ALIGN_RNA_FASTP {

    take:
        reads_ch
        star_index

    main:
        fp  = FASTP(reads_ch)
        aln = STAR_RNA_ALIGN(fp.reads, star_index)

    emit:
        bam         = aln.bam
        sj          = aln.sj
        gene_counts = aln.gene_counts
        star_log    = aln.log

        fastp_json  = fp.json
        fastp_html  = fp.html
        cleaned_reads = fp.reads
}