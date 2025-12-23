nextflow.enable.dsl=2

include { FASTP } from '../modules/fastp/main.nf'

workflow QC_FASTP {

    take:
        reads_ch   // tuple(sample_id, fastq)

    main:
        FASTP(reads_ch)

    emit:
        trimmed_reads = FASTP.out.trimmed_reads
}
