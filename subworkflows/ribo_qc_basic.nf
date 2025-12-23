nextflow.enable.dsl=2

include { RPF_LENGTH_QC } from '../modules/rpf_length/main.nf'

workflow RIBO_QC_BASIC {

    take:
        aligned_ch   // tuple(sample_id, bam, bai)

    main:
        // Explicitly preserve BAM + BAI
        qc_input = aligned_ch.map { sample_id, bam, bai ->
            tuple(sample_id, bam, bai)
        }

        RPF_LENGTH_QC(qc_input)

    emit:
        RPF_LENGTH_QC.out
}