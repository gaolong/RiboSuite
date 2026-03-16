nextflow.enable.dsl = 2

process STAR_RNA_ALIGN {

    tag "${meta.sample_id}"

    conda "bioconda::star=2.7.11b"

    publishDir "${params.outdir}/rna/align/star_rna",
        mode: 'copy',
        saveAs: { file -> "${meta.sample_id}/${file}" }

    input:
    tuple val(meta), path(reads), path(star_index)

    output:
    tuple val(meta), path("${meta.sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.sample_id}.SJ.out.tab"),                    emit: sj
    tuple val(meta), path("${meta.sample_id}.ReadsPerGene.out.tab"),          emit: gene_counts
    tuple val(meta), path("${meta.sample_id}.Log.final.out"),                 emit: log

    script:
    def read_list = (reads instanceof List) ? reads : [reads]
    def read_args = read_list.join(' ')

    """
    STAR \
      --genomeDir ${star_index} \
      --readFilesIn ${read_args} \
      --readFilesCommand zcat \
      --runThreadN ${task.cpus} \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix ${meta.sample_id}.
    """
}