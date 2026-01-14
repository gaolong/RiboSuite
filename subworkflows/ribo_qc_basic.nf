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
         */
        RPF_LENGTH_QC(aligned_ch)

        /*
         * 2) P-site offset
         * Emits: (sample_id, bam, bai, offsets)
         */
        PSITE_OFFSET(aligned_ch, gtf)

        /*
         * 3) Frame periodicity QC
         * Consumes: (sample_id, bam, bai, offsets)
         */
        FRAME_PERIODICITY_QC(
            PSITE_OFFSET.out,
            gtf
        )

        /*
         * 4) Metagene QC
         * MUST receive offsets
         */
        METAGENE_QC(
            PSITE_OFFSET.out.map { sample_id, bam, bai, offsets ->
                tuple(sample_id, bam, bai, offsets)
            },
            gtf
        )

    emit:
        /*
         * QC tables
         */
        rpf_length_qc         = RPF_LENGTH_QC.out
        psite_offset_qc       = PSITE_OFFSET.out

        /*
         * Frame periodicity (split outputs)
         */
        frame_periodicity_tsv = FRAME_PERIODICITY_QC.out.periodicity_tsv
        frame_periodicity_png = FRAME_PERIODICITY_QC.out.periodicity_png

        /*
         * Metagene QC
         */
        metagene_qc           = METAGENE_QC.out
}
