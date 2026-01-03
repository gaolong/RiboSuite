// DEPRECATED: RNA-seq trimming is now handled by FASTP in rna_preprocess.nf
// This subworkflow should not be used for RNA-seq

include { CUTADAPT } from '../modules/cutadapt/main.nf'

workflow RNA_TRIM_READS {

    take:
    reads_ch

    main:
    reads_ch | CUTADAPT

    emit:
    CUTADAPT.out  as trimmed_reads
}
