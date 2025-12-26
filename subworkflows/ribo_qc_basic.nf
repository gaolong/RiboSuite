nextflow.enable.dsl=2

include { RPF_LENGTH_QC }        from '../modules/rpf_length/main.nf'
include { PSITE_OFFSET }         from '../modules/psite_offset/main.nf'
include { FRAME_PERIODICITY_QC } from '../modules/frame_periodicity/main.nf'
include { METAGENE_QC }          from '../modules/metagene/main.nf'


workflow RIBO_QC_BASIC {

    take:
        aligned_ch   // tuple(sample_id, bam, bai)
        gtf          // path to GTF file

    main:

        /*
         * 1) Read-length QC
         * Consumes (sample_id, bam, bai)
         * bai is staged but not used
         */
        RPF_LENGTH_QC(aligned_ch)

        /*
         * 2) P-site offset
         * Consumes (sample_id, bam, bai) + gtf
         * Emits (sample_id, bam, bai, psite_offsets.tsv)
         */
        PSITE_OFFSET(aligned_ch, gtf)

        /*
         * 3) Frame periodicity QC
         * Consumes (sample_id, bam, bai, psite_offsets.tsv) + gtf
         */
        FRAME_PERIODICITY_QC(PSITE_OFFSET.out, gtf)

        /*
         * 4) Metagene QC
         * Uses P-site–corrected BAM
         */
        METAGENE_QC(
            PSITE_OFFSET.out.map { sample_id, bam, bai, offsets ->
                tuple(sample_id, bam, bai)
            },
            gtf
        )

    emit:
        RPF_LENGTH_QC.out
        PSITE_OFFSET.out
        FRAME_PERIODICITY_QC.out
        METAGENE_QC.out
}
  
