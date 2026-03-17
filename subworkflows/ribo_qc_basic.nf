nextflow.enable.dsl = 2

include { RPF_LENGTH_QC }        from '../modules/rpf_length/main.nf'
include { PSITE_OFFSET }         from '../modules/psite_offset/main.nf'
include { FRAME_PERIODICITY_QC } from '../modules/frame_periodicity/main.nf'
include { METAGENE_QC }          from '../modules/metagene/main.nf'

workflow RIBO_QC_BASIC {

    take:
        // tuple(meta, bam, bai, gtf)
        aligned_ch

    main:

        /*
         * 1) Read-length QC
         * RPF_LENGTH_QC likely still expects: (sample_id, bam, bai)
         */
        rpf_len_in = aligned_ch.map { meta, bam, bai, gtf ->
            tuple(meta.sample_id.toString(), bam, bai)
        }

        RPF_LENGTH_QC(rpf_len_in)

        /*
         * 2) P-site offset
         * New expected input:
         *   tuple(meta, bam, bai, gtf)
         * New expected output:
         *   tuple(meta, bam, bai, offsets)
         */
        psite_in = aligned_ch.map { meta, bam, bai, gtf ->
            tuple(meta, bam, bai, gtf)
        }

        PSITE_OFFSET(psite_in)

        /*
         * 3) Frame periodicity QC
         * New expected input:
         *   tuple(meta, bam, bai, offsets, gtf)
         */
        frame_in = PSITE_OFFSET.out
            .join(
                aligned_ch.map { meta, bam, bai, gtf -> tuple(meta, gtf) },
                by: 0
            )
            .map { meta, bam, bai, offsets, gtf ->
                tuple(meta, bam, bai, offsets, gtf)
            }

        FRAME_PERIODICITY_QC(frame_in)

        /*
         * 4) Metagene QC
         * New expected input:
         *   tuple(meta, bam, bai, offsets, gtf)
         */
        metagene_in = PSITE_OFFSET.out
            .join(
                aligned_ch.map { meta, bam, bai, gtf -> tuple(meta, gtf) },
                by: 0
            )
            .map { meta, bam, bai, offsets, gtf ->
                tuple(meta, bam, bai, offsets, gtf)
            }

        METAGENE_QC(metagene_in)

    emit:
        /*
         * QC tables
         */
        rpf_length_qc   = RPF_LENGTH_QC.out
        psite_offset_qc = PSITE_OFFSET.out

        /*
         * Frame periodicity (split outputs)
         */
        frame_periodicity_tsv = FRAME_PERIODICITY_QC.out.periodicity_tsv
        frame_periodicity_png = FRAME_PERIODICITY_QC.out.periodicity_png

        /*
         * Metagene QC
         */
        metagene_qc = METAGENE_QC.out
}