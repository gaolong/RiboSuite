nextflow.enable.dsl = 2

include { STAR_RIBO_ALIGN }           from '../modules/star_ribo/main.nf'
include { SAMTOOLS_INDEX }            from '../modules/samtools/index/main.nf'
include { SAMTOOLS_FILTER_BY_QNAME }  from '../modules/samtools/filter_by_qname/main.nf'
include { UMI_DEDUP }                 from '../modules/umi_dedup/main.nf'

workflow ALIGN_RIBO_DEDUP {

    take:
        reads_ch        // (sid, clean_fastq, has_umi, sjdb_or_null_or_NO_FILE)
        star_index

    main:

        /*
         * Sentinel used by STAR_RIBO_ALIGN to mean "no sjdb"
         * Make sure this file exists: assets/NO_FILE
         */
        def NO_FILE = file("${projectDir}/assets/NO_FILE")

        /*
         * 0) Side metadata
         */
        meta_side = reads_ch.map { sid, fq, has_umi, sjdb ->
            tuple(sid, has_umi, sjdb)
        }

        /*
         * Normalize sjdb: null -> NO_FILE
         */
        filtered = reads_ch.map { sid, fq, has_umi, sjdb ->
            def sjdb_path = (sjdb != null) ? sjdb : NO_FILE
            tuple(sid, fq, has_umi, sjdb_path)
        }
        // (sid, clean_fastq, has_umi, sjdb_path)

        /*
         * 1) STAR alignment
         */
        STAR_RIBO_ALIGN(
            filtered.map { sid, fq, has_umi, sjdb_path -> tuple(sid, fq) },
            file(star_index),
            filtered.map { sid, fq, has_umi, sjdb_path -> sjdb_path }
        )

        aligned_genome = STAR_RIBO_ALIGN.out.genome_bam
            .join(meta_side.map { sid, has_umi, sjdb -> tuple(sid, has_umi) }, by: 0)
        // (sid, genome_bam, has_umi)

        aligned_tx = STAR_RIBO_ALIGN.out.tx_bam

        /*
         * 2) Index
         */
        indexed0 = SAMTOOLS_INDEX(
            aligned_genome.map { sid, genome_bam, has_umi -> tuple(sid, genome_bam) }
        )
        indexed = indexed0.join(meta_side.map { sid, has_umi, sjdb -> tuple(sid, has_umi) }, by: 0)
        // (sid, genome_bam, bai, has_umi)

        /*
         * 3) UMI dedup
         */
        with_umi    = indexed.filter { sid, bam, bai, has_umi -> has_umi }
        without_umi = indexed.filter { sid, bam, bai, has_umi -> !has_umi }

        deduped = UMI_DEDUP(
            with_umi.map { sid, bam, bai, has_umi -> tuple(sid, bam, bai) }
        )
        // (sid, dedup_bam, dedup_bai, dedup_log, survivors_qnames)

        deduped_genome = deduped.map { sid, dedup_bam, dedup_bai, log, qnames ->
            tuple(sid, dedup_bam, dedup_bai)
        }

        passed_genome = without_umi.map { sid, bam, bai, has_umi ->
            tuple(sid, bam, bai)
        }

        genome_bam_out = deduped_genome.mix(passed_genome)

        /*
         * 4) tx filtering (UMI samples only)
         */
        tx_with_umi = aligned_tx
            .join(deduped, by: 0)
            .map { sid, tx_bam, dedup_bam, dedup_bai, dedup_log, qnames ->
                tuple(sid, tx_bam, qnames)
            }

        tx_filtered = SAMTOOLS_FILTER_BY_QNAME(tx_with_umi)

    emit:
        genome_bam = genome_bam_out
        tx_bam     = tx_filtered
}