#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * RNA-seq quantification workflow
 * - raw QC
 * - trimming
 * - STAR alignment
 * - gene-level quantification
 */

include { RNA_QC_RAW }      from '../subworkflows/rna_qc_raw.nf'
include { RNA_TRIM_READS }  from '../subworkflows/rna_trim_reads.nf'
include { ALIGN_RNA }       from '../subworkflows/align_rna.nf'
include { FEATURECOUNTS }   from '../modules/subread/featurecounts/main.nf'


workflow RNA_QUANT {

    /*
     * Input:
     *   samplesheet with sample_id + fastq_1 (+ fastq_2 optional)
     */
    main:

    /*
     * Parse samplesheet
     * Expected output:
     *   tuple(meta, reads)
     *     meta.sample_id
     *     reads = [fq] or [fq1, fq2]
     */
    reads_ch = Channel
        .fromPath(params.samplesheet)
        | parseSamplesheet


    /*
     * 1. Raw QC (FastQC only, non-destructive)
     */
    qc_out = RNA_QC_RAW(reads_ch)


    /*
     * 2. Adapter trimming
     */
    trimmed_ch = RNA_TRIM_READS(qc_out.reads_raw)


    /*
     * 3. Alignment (STAR → sort → index)
     */
    aligned_ch = ALIGN_RNA(trimmed_ch.trimmed_reads)


    /*
     * 4. Gene-level quantification
     */
    FEATURECOUNTS(
        aligned_ch.aligned_bam,
        params.gtf
    )
}
