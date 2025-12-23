nextflow.enable.dsl=2

include { RPF_LENGTH_QC } from '../modules/rpf_length/main.nf'

workflow RIBO_QC_BASIC {

    take:
        bam_ch   // tuple(sample_id, sorted.bam)

    main:
        rpf_len = RPF_LENGTH_QC(bam_ch)

    emit:
        rpf_len
}