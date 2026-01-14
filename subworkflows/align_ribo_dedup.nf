nextflow.enable.dsl = 2

include { BOWTIE2_FILTER }            from '../modules/bowtie2/main.nf'
include { STAR_RIBO_ALIGN }          from '../modules/star_ribo/main.nf'
include { SAMTOOLS_INDEX }           from '../modules/samtools/index/main.nf'
include { SAMTOOLS_FILTER_BY_QNAME } from '../modules/samtools/filter_by_qname/main.nf'
include { UMI_DEDUP }                from '../modules/umi_dedup/main.nf'


workflow ALIGN_RIBO_DEDUP {

    take:
        reads_ch        // (sid, clean_fastq, has_umi)
        contam_index
        star_index

    main:

        /*
         * ------------------------------------------------------------
         * 0. Side metadata channel
         * ------------------------------------------------------------
         */
        meta_has_umi = reads_ch.map { sid, fq, has_umi ->
            tuple(sid, has_umi)
        }

        /*
         * ------------------------------------------------------------
         * 1. Contaminant filtering
         *    Input : (sid, fastq)
         *    Output: (sid, clean_fastq)
         * ------------------------------------------------------------
         */
        filtered_reads = BOWTIE2_FILTER(
            reads_ch.map { sid, fq, has_umi -> tuple(sid, fq) },
            file(contam_index).parent,
            file(contam_index).getName()
        ).clean_reads

        filtered = filtered_reads.join(meta_has_umi)
        // (sid, clean_fastq, has_umi)

        /*
         * ------------------------------------------------------------
         * 2. STAR genome / transcriptome alignment
         * ------------------------------------------------------------
         */
        STAR_RIBO_ALIGN(
            filtered.map { sid, fq, has_umi -> tuple(sid, fq) },
            file(star_index)
        )

        /*
         * Explicit STAR outputs
         */
        aligned_genome = STAR_RIBO_ALIGN.out.genome_bam
            .join(meta_has_umi)
        // (sid, genome_bam, has_umi)

        aligned_tx = STAR_RIBO_ALIGN.out.tx_bam
        // (sid, tx_bam)   [optional, may be empty]

        /*
         * ------------------------------------------------------------
         * 3. Index genome BAM
         * ------------------------------------------------------------
         */
        indexed0 = SAMTOOLS_INDEX(
            aligned_genome.map { sid, genome_bam, has_umi ->
                tuple(sid, genome_bam)
            }
        )
        // (sid, genome_bam, bai)

        indexed = indexed0.join(meta_has_umi)
        // (sid, genome_bam, bai, has_umi)

        /*
         * ------------------------------------------------------------
         * 4. Genome-space UMI deduplication
         * ------------------------------------------------------------
         */
        with_umi    = indexed.filter { sid, bam, bai, has_umi -> has_umi }
        without_umi = indexed.filter { sid, bam, bai, has_umi -> !has_umi }

        deduped = UMI_DEDUP(
            with_umi.map { sid, bam, bai, has_umi ->
                tuple(sid, bam, bai)
            }
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
         * ------------------------------------------------------------
         * 5. Best-effort transcriptome BAM filtering
         *    (UMI samples only, quantMode enabled only)
         * ------------------------------------------------------------
         */
        tx_with_umi = aligned_tx
            .join(deduped)
            .map { sid,
                   tx_bam,
                   dedup_bam,
                   dedup_bai,
                   dedup_log,
                   qnames ->
                tuple(sid, tx_bam, qnames)
            }

        tx_filtered = SAMTOOLS_FILTER_BY_QNAME(tx_with_umi)
        // (sid, tx_filtered_bam, tx_filtered_bam.bai)

    emit:
        /*
         * Canonical genome BAM for QC & downstream analysis
         */
        genome_bam = genome_bam_out

        /*
         * Optional transcriptome BAM
         * (only when --star_quantmode true AND has_umi)
         */
        tx_bam     = tx_filtered
}
