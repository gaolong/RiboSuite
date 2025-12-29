nextflow.enable.dsl=2

include { PREPROCESS_READS } from '../subworkflows/preprocess_reads.nf'
include { ALIGN_RIBO_DEDUP } from '../subworkflows/align_ribo_dedup.nf'
include { RIBO_QC_BASIC }    from '../subworkflows/ribo_qc_basic.nf'

workflow TEST_PREPROCESS_ALIGN_QC {

    main:
        /*
         * Load samples
         */
        samples = Channel
            .fromPath(params.samples)
            .splitCsv(header:true, sep:'\t')
            .map { row ->
                def has_umi = row.bc_pattern && row.bc_pattern.trim()
                tuple(
                    row.sample_id,
                    file(row.fastq),
                    row.adapter,
                    row.bc_pattern,
                    has_umi
                )
            }

        /*
         * 1. Preprocess (cutadapt → umi_extract → fastqc)
         * Emits: (sample_id, clean_fastq, has_umi)
         */
        preprocessed = PREPROCESS_READS(samples)

        /*
         * 2. Alignment + optional UMI dedup
         */
        bam = ALIGN_RIBO_DEDUP(
            preprocessed,
            params.contam_index,
            params.star_index
        )

        /*
         * 3. Ribo QC
         */
        qc = RIBO_QC_BASIC(
            bam,
            file(params.gtf)
        )

    emit:
        rpf_length_qc        = qc.rpf_length_qc
        psite_offset_qc      = qc.psite_offset_qc
        frame_periodicity_qc = qc.frame_periodicity_qc
        metagene_qc          = qc.metagene_qc
}
