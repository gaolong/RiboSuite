nextflow.enable.dsl=2

include { QC_FASTP } from '../subworkflows/qc_fastp.nf'

workflow RiboSuite {
    take:
        reads = params.reads

    main:
        trimmed = QC_FASTP(reads)

    emit:
        trimmed
}