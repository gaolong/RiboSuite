include { CUTADAPT } from '../modules/cutadapt/main.nf'

workflow RNA_TRIM_READS {

    take:
    reads_ch

    main:
    reads_ch | CUTADAPT

    emit:
    CUTADAPT.out  as trimmed_reads
}
