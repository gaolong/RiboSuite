nextflow.enable.dsl = 2

process STAR_RNA_ALIGN {

    tag "${meta.sample_id}"

    conda "bioconda::star=2.7.11b"

    publishDir "${params.outdir}/rna/align/star_ribo",
               mode: 'copy',
               saveAs: { file -> "${sample_id}/${file}" }

    input:
    tuple val(meta), path(read1), path(read2)
    path star_index

    output:
    tuple val(meta),
          path("${meta.sample_id}.Aligned.sortedByCoord.out.bam"),
          emit: bam

    path "${meta.sample_id}.Log.final.out", emit: log
    path "${meta.sample_id}.ReadsPerGene.out.tab", emit: gene_counts
    path "${meta.sample_id}.SJ.out.tab", emit: sj

    script:
    """
    STAR \
      --genomeDir ${star_index} \
      --readFilesIn ${read1} ${read2} \
      --readFilesCommand zcat \
      --runThreadN ${task.cpus} \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix ${meta.sample_id}.
    """
}