nextflow.enable.dsl=2

include { RPF_LENGTH_QC } from '../modules/rpf_length/main.nf'
include { PSITE_OFFSET }   from '../modules/psite_offset/main.nf'

workflow RIBO_QC_BASIC {

    take:
        aligned_ch   // tuple(sample_id, bam)
        gtf

    main:
        // preserve tuple shape explicitly
        qc_input = aligned_ch.map { sample_id, bam ->
            tuple(sample_id, bam)
        }

        RPF_LENGTH_QC(qc_input)
        PSITE_OFFSET(qc_input, gtf)

    emit:
        RPF_LENGTH_QC.out
        PSITE_OFFSET.out
}
