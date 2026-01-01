include { STAR } from '../modules/star/main.nf'
include { SAMTOOLS_SORT } from '../modules/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/samtools/index/main.nf'

workflow ALIGN_RNA {

    take:
    reads_ch

    main:
    reads_ch | STAR | SAMTOOLS_SORT | SAMTOOLS_INDEX

    emit:
    SAMTOOLS_INDEX.out as aligned_bam
}
