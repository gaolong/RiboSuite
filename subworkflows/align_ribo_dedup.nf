nextflow.enable.dsl = 2

include { STAR_RIBO_ALIGN }           from '../modules/star_ribo/main.nf'
include { SAMTOOLS_INDEX }            from '../modules/samtools/index/main.nf'
include { SAMTOOLS_FILTER_BY_QNAME }  from '../modules/samtools/filter_by_qname/main.nf'
include { UMI_DEDUP }                 from '../modules/umi_dedup/main.nf'

workflow ALIGN_RIBO_DEDUP {

    take:
        // tuple(meta, clean_fastq, has_umi, sjdb_or_null_or_NO_FILE, star_index)
        reads_ch

    main:

        /*
         * Sentinel used by STAR_RIBO_ALIGN to mean "no sjdb"
         */
        def NO_FILE = file("${projectDir}/assets/NO_FILE")

        /*
         * 0) Side metadata
         */
        meta_side = reads_ch.map { meta, fq, has_umi, sjdb, star_index ->
            tuple(meta, has_umi, sjdb)
        }

        /*
         * Normalize sjdb: null -> NO_FILE
         */
        filtered = reads_ch.map { meta, fq, has_umi, sjdb, star_index ->
            def sjdb_path = (sjdb != null) ? sjdb : NO_FILE
            tuple(meta, fq, has_umi, sjdb_path, star_index)
        }
        // (meta, clean_fastq, has_umi, sjdb_path, star_index)

        /*
         * 1) STAR alignment
         */
        STAR_RIBO_ALIGN(filtered)

        aligned_genome = STAR_RIBO_ALIGN.out.genome_bam
            .join(
                meta_side.map { meta, has_umi, sjdb -> tuple(meta, has_umi) },
                by: 0
            )
        // (meta, genome_bam, has_umi)

        aligned_tx = STAR_RIBO_ALIGN.out.tx_bam

        /*
         * 2) Index
         */
        indexed0 = SAMTOOLS_INDEX(
            aligned_genome.map { meta, genome_bam, has_umi ->
                tuple(meta.sample_id.toString(), genome_bam)
            }
        )

        indexed = indexed0
            .map { sid, bam, bai -> tuple(sid, bam, bai) }
            .join(
                meta_side.map { meta, has_umi, sjdb ->
                    tuple(meta.sample_id.toString(), meta, has_umi)
                },
                by: 0
            )
            .map { sid, bam, bai, meta, has_umi ->
                tuple(meta, bam, bai, has_umi)
            }
        // (meta, genome_bam, bai, has_umi)

        /*
         * 3) UMI dedup
         */
        with_umi    = indexed.filter  { meta, bam, bai, has_umi -> has_umi }
        without_umi = indexed.filter  { meta, bam, bai, has_umi -> !has_umi }

        deduped = UMI_DEDUP(
            with_umi.map { meta, bam, bai, has_umi ->
                tuple(meta.sample_id.toString(), bam, bai)
            }
        )
        // (sid, dedup_bam, dedup_bai, dedup_log, survivors_qnames)

        deduped_genome = deduped
            .join(
                with_umi.map { meta, bam, bai, has_umi ->
                    tuple(meta.sample_id.toString(), meta)
                },
                by: 0
            )
            .map { sid, dedup_bam, dedup_bai, log, qnames, meta ->
                tuple(meta, dedup_bam, dedup_bai)
            }

        passed_genome = without_umi.map { meta, bam, bai, has_umi ->
            tuple(meta, bam, bai)
        }

        genome_bam_out = deduped_genome.mix(passed_genome)

        /*
         * 4) tx filtering (UMI samples only)
         */
        tx_with_umi = aligned_tx
            .map { meta, tx_bam -> tuple(meta.sample_id.toString(), meta, tx_bam) }
            .join(deduped, by: 0)
            .map { sid, meta, tx_bam, dedup_bam, dedup_bai, dedup_log, qnames ->
                tuple(sid, tx_bam, qnames)
            }

        tx_filtered0 = SAMTOOLS_FILTER_BY_QNAME(tx_with_umi)

        tx_filtered = tx_filtered0
            .join(
                with_umi.map { meta, bam, bai, has_umi ->
                    tuple(meta.sample_id.toString(), meta)
                },
                by: 0
            )
            .map { sid, tx_bam, meta ->
                tuple(meta, tx_bam)
            }

    emit:
        genome_bam = genome_bam_out
        tx_bam     = tx_filtered
}