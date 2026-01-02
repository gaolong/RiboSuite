include { STAR_RNA_ALIGN }  from '../modules/star_rna/main.nf'
include { SAMTOOLS_INDEX }  from '../modules/samtools/index/main.nf'

workflow ALIGN_RNA {

    take:
    reads_ch

    main:
    /*
     * STAR requires TWO inputs:
     *  1) reads_ch
     *  2) star index
     */
    STAR_RNA_ALIGN(reads_ch, params.star_index)

    SAMTOOLS_INDEX(STAR_RNA_ALIGN.out.bam)

    emit:
    SAMTOOLS_INDEX.out
    STAR_RNA_ALIGN.out.gene_counts
    STAR_RNA_ALIGN.out.log
}
