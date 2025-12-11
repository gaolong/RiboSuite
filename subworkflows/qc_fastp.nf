nextflow.enable.dsl=2

include { FASTP } from '../modules/fastp/main.nf'

workflow QC_FASTP {

    take:
        reads_ch

    main:
        trimmed_ch = FASTP(reads_ch)

    emit:
        trimmed_ch
}