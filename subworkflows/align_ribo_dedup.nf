nextflow.enable.dsl=2

include { BOWTIE2_FILTER }  from '../modules/bowtie2/main.nf'
include { STAR_ALIGN }      from '../modules/star/main.nf'
include { SAMTOOLS_INDEX }  from '../modules/samtools/index/main.nf'
include { UMI_DEDUP }       from '../modules/umi_dedup/main.nf'

workflow ALIGN_RIBO_DEDUP {

    take:
        reads_ch        // (sample_id, clean_fastq, has_umi)
        contam_index
        star_index

    main:
        /*
         * Keep metadata as a side channel keyed by sample_id
         */
        meta_has_umi = reads_ch.map { sid, fq, has_umi -> tuple(sid, has_umi) }

        /*
         * 1. Contaminant filtering
         * Feed bowtie only (sid, fastq)
         */
        filtered0 = BOWTIE2_FILTER(
            reads_ch.map { sid, fq, has_umi -> tuple(sid, fq) },
            file(contam_index).parent,
            file(contam_index).getName()
        ).clean_reads   // expected: (sid, clean_fastq)

        filtered = filtered0.join(meta_has_umi) // (sid, clean_fastq, has_umi)

        /*
         * 2. STAR alignment
         * Feed STAR only (sid, fastq)
         */
        aligned0 = STAR_ALIGN(
            filtered.map { sid, fq, has_umi -> tuple(sid, fq) },
            file(star_index)
        )  // expected: (sid, bam)

        aligned = aligned0.join(meta_has_umi) // (sid, bam, has_umi)

        /*
         * 3. BAM indexing
         */
        indexed0 = SAMTOOLS_INDEX(
            aligned.map { sid, bam, has_umi -> tuple(sid, bam) }
        ) // expected: (sid, bam, bai)

        indexed = indexed0.join(meta_has_umi) // (sid, bam, bai, has_umi)

        /*
         * 4. Conditional UMI dedup
         */
        with_umi    = indexed.filter { sid, bam, bai, has_umi -> has_umi }
        without_umi = indexed.filter { sid, bam, bai, has_umi -> !has_umi }

        deduped = UMI_DEDUP(
            with_umi.map { sid, bam, bai, has_umi -> tuple(sid, bam, bai) }
        )
        .map { sid, bam, bai, log ->
            tuple(sid, bam, bai)
        }

        passed = without_umi.map { sid, bam, bai, has_umi ->
            tuple(sid, bam, bai)
        }

        bam_out = deduped.mix(passed)

    emit:
        bam = bam_out
}
