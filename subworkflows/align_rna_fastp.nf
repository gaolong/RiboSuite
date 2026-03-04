nextflow.enable.dsl = 2

include { FASTP }          from '../modules/fastp/main.nf'
include { STAR_RNA_ALIGN } from '../modules/star_rna/main.nf'

workflow ALIGN_RNA_FASTP {

    take:
        reads_ch
        star_index

    main:
        fp  = FASTP(reads_ch)
        aln = STAR_RNA_ALIGN(fp.reads, star_index)

        // ------------------------------------------------------------
        // Normalize outputs to Nextflow tuples (avoid java.util.ArrayList)
        // ------------------------------------------------------------
        bam_ch = aln.bam.map { row -> tuple(row[0], row[1]) }

        sj_ch  = aln.sj.map  { row -> tuple(row[0], row[1]) }
        gc_ch  = aln.gene_counts.map { row -> tuple(row[0], row[1]) }
        log_ch = aln.log.map { row -> tuple(row[0], row[1]) }

    emit:
        bam         = bam_ch
        sj          = sj_ch
        gene_counts = gc_ch
        star_log    = log_ch

        fastp_json  = fp.json
        fastp_html  = fp.html
        cleaned_reads = fp.reads
}