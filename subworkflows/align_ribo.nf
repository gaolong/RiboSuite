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
        filtered = BOWTIE2_FILTER(reads_ch, contam_index)
        aligned  = STAR_ALIGN(filtered.map { it[0], it[1] }, star_index)
        sorted   = SAMTOOLS_SORT_INDEX(aligned)

    emit:
        sorted
}