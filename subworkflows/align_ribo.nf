nextflow.enable.dsl=2

include { BOWTIE2_FILTER }      from '../modules/bowtie2/main.nf'
include { STAR_ALIGN }          from '../modules/star/main.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools/main.nf'

workflow ALIGN_RIBO {

    take:
        reads_ch        // tuple(sample_id, trimmed_fastq)
        contam_index    // val (string path prefix)
        star_index      // val (STAR index dir)

    main:
        // Contamination filtering (rRNA/tRNA/etc)
        filtered = BOWTIE2_FILTER(
            reads_ch,
            contam_index
        )

        // STAR alignment on clean reads
        aligned = STAR_ALIGN(
            filtered.clean_reads,
            star_index
        )

        // Sort & index BAM
        sorted = SAMTOOLS_SORT_INDEX(aligned)

    emit:
        sorted
}
