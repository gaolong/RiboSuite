nextflow.enable.dsl=2

include { BOWTIE2_FILTER } from '../modules/bowtie2/main.nf'
include { STAR_ALIGN }     from '../modules/star/main.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools/main.nf'

workflow ALIGN_RIBO {

    take:
        reads_ch
        contam_index
        star_index

    main:
        // Bowtie2 contamination filtering
        filtered = BOWTIE2_FILTER(
            reads_ch.map { sample_id, fastq, json, html -> tuple(sample_id, fastq) },
            contam_index
        )

        // Use ONLY the clean reads for STAR
        aligned = STAR_ALIGN(filtered.clean_reads, star_index)

        // Sort & index BAM
        sorted = SAMTOOLS_SORT_INDEX(aligned)

    emit:
        sorted
}
