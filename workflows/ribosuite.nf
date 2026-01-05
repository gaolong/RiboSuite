nextflow.enable.dsl = 2

/*
 * Subworkflows
 */
include { PREPROCESS_READS }   from '../subworkflows/preprocess_reads.nf'
include { ALIGN_RIBO_DEDUP }   from '../subworkflows/align_ribo_dedup.nf'
include { RIBO_QC_BASIC }      from '../subworkflows/ribo_qc_basic.nf'
include { CDS_QUANT }          from '../modules/cds_quant/main.nf'


/*
 * ------------------------------------------------------------
 * Core pipeline workflow (logic only)
 * ------------------------------------------------------------
 */
workflow RiboSuite {

    take:
        reads_ch    // (sample_id, fastq, adapter, bc_pattern, has_umi)
        gtf         // annotation GTF

    main:

        /*
         * 1) Preprocess reads
         * Emits: (sample_id, clean_fastq, has_umi)
         */
        processed = PREPROCESS_READS(reads_ch)

        /*
         * 2) Alignment + optional UMI dedup
         * Emits: (sample_id, bam, bai)
         */
        aligned = ALIGN_RIBO_DEDUP(
            processed,
            params.contam_index,
            params.star_index
        )

        /*
         * 3) Ribo-seq QC
         * Includes PSITE_OFFSET internally
         */
        qc = RIBO_QC_BASIC(
            aligned,
            gtf
        )

        /*
         * 4) CDS quantification
         * MUST use P-site offsets
         */

        // Extract (sample_id, bam, bai, offsets)
        psite_for_quant = qc.psite_offset_qc.map { sid, bam, bai, offsets ->
            tuple(sid, bam, bai, offsets)
        }

        cds_quant = CDS_QUANT(
            // (sample_id, bam, bai)
            psite_for_quant.map { sid, bam, bai, offsets ->
                tuple(sid, bam, bai)
            },
            // gtf
            gtf,
            // offsets
            psite_for_quant.map { sid, bam, bai, offsets ->
                offsets
            }
        )

    emit:
        aligned
        rpf_length_qc        = qc.rpf_length_qc
        psite_offset_qc      = qc.psite_offset_qc
        frame_periodicity_qc = qc.frame_periodicity_qc
        metagene_qc          = qc.metagene_qc
        cds_quant
}


/*
 * ------------------------------------------------------------
 * Default entry workflow (CLI-facing wrapper)
 * ------------------------------------------------------------
 */
workflow {

    /*
     * Build reads channel from samples.tsv
     */
    reads_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->

            def has_umi = row.bc_pattern && row.bc_pattern != 'NA'

            tuple(
                row.sample_id,
                file(row.fastq),
                row.adapter,
                row.bc_pattern,
                has_umi
            )
        }

    /*
     * Launch core workflow
     */
    RiboSuite(
        reads_ch,
        file(params.gtf)
    )
}
