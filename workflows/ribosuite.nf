nextflow.enable.dsl = 2

/*
 * ------------------------------------------------------------
 * Subworkflows / modules
 * ------------------------------------------------------------
 */
include { PREPROCESS_READS }   from '../subworkflows/preprocess_reads.nf'
include { ALIGN_RIBO_DEDUP }   from '../subworkflows/align_ribo_dedup.nf'
include { RIBO_QC_BASIC }      from '../subworkflows/ribo_qc_basic.nf'
include { CDS_QUANT }          from '../modules/cds_quant/main.nf'
include { PSITE_TRACK }        from '../modules/psite_track/main.nf'
include { TRANSLATED_ORFS }    from '../modules/translated_orfs/main.nf'


/*
 * ------------------------------------------------------------
 * Core pipeline workflow (logic only)
 * ------------------------------------------------------------
 */
workflow RiboSuite {

    take:
        // (sample_id, fastq, adapter_5, adapter_3, bc_pattern, has_umi)
        reads_ch
        gtf

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
         *
         * qc.psite_offset_qc emits: (sample_id, bam, bai, offsets)
         */
        qc = RIBO_QC_BASIC(
            aligned.genome_bam,
            gtf
        )

        /*
         * 4) CDS quantification (requires P-site offsets)
         */
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

        /*
         * 5) Optional translated ORF calling (uses bam + offsets + gtf)
         */
        if (params.enable_translated_orfs) {
            translated_orfs = TRANSLATED_ORFS(
                // (sample_id, bam, bai, offsets)
                qc.psite_offset_qc.map { sid, bam, bai, offsets ->
                    tuple(sid, bam, bai, offsets)
                },
                // gtf
                gtf
            )
        }

        /*
         * 6) Optional P-site bigWig tracks
         */
        if (params.enable_psite_track) {

            psite_track = PSITE_TRACK(
                // (sample_id, bam)
                aligned.genome_bam.map { sid, bam, bai ->
                    tuple(sid, bam)
                },
                // offsets
                qc.psite_offset_qc.map { sid, bam, bai, offsets ->
                    offsets
                },
                // genome sizes
                file(params.genome_sizes)
            )
        }

    emit:
        genome_bam              = aligned.genome_bam

        rpf_length_qc           = qc.rpf_length_qc
        psite_offset_qc         = qc.psite_offset_qc

        frame_periodicity_tsv   = qc.frame_periodicity_tsv
        frame_periodicity_png   = qc.frame_periodicity_png

        metagene_qc             = qc.metagene_qc

        cds_quant

        
        translated_orfs_tsv = params.enable_translated_orfs \
            ? translated_orfs.tsv \
            : Channel.empty()

        psite_tracks = params.enable_psite_track \
            ? psite_track.psite_tracks_by_len \
            : Channel.empty()
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

            /*
             * Normalize adapters
             */
            def adapter_5 = (row.adapter_5 && row.adapter_5 != '0' && row.adapter_5 != 'NA')
                ? row.adapter_5
                : null

            def adapter_3 = (row.adapter_3 && row.adapter_3 != '0' && row.adapter_3 != 'NA')
                ? row.adapter_3
                : null

            /*
             * Normalize UMI lengths
             */
            int umi5 = row.umi_5?.trim() ? row.umi_5.toInteger() : 0
            int umi3 = row.umi_3?.trim() ? row.umi_3.toInteger() : 0

            /*
             * Generate UMI bc_pattern
             */
            def bc_pattern = null
            if (umi5 > 0 || umi3 > 0) {
                def parts = []
                if (umi5 > 0) parts << "(?P<umi_1>.{${umi5}})"
                parts << ".+"
                if (umi3 > 0) parts << "(?P<umi_2>.{${umi3}})"
                bc_pattern = "^${parts.join('')}\$"
            }

            def has_umi = bc_pattern != null

            /*
             * Optional sanity warning
             */
            if (has_umi && adapter_5 == null && adapter_3 == null) {
                log.warn "Sample ${row.sample_id}: UMI specified but no adapters provided"
            }

            tuple(
                row.sample_id,
                file(row.fastq),
                adapter_5,
                adapter_3,
                bc_pattern,
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