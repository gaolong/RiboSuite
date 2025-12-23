nextflow.enable.dsl=2

include { QC_FASTP } from '../subworkflows/qc_fastp.nf'
include { ALIGN_RIBO } from '../subworkflows/align_ribo.nf'
include { RIBO_QC_BASIC } from '../subworkflows/ribo_qc_basic.nf'

workflow RiboSuite {

    take:
        reads   // just a name, no assignment

    main:
        trimmed = QC_FASTP(reads)

        aligned = ALIGN_RIBO(
            trimmed,
            params.contam_index,
            params.star_index
        )

        qc_basic = RIBO_QC_BASIC(aligned)

    emit:
        qc_basic
}