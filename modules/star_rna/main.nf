nextflow.enable.dsl = 2

process STAR_RNA_ALIGN {

    tag "${meta.sample_id}"

    conda "bioconda::star=2.7.11b"

    publishDir "${params.outdir}/rna/align/star_rna",
           mode: 'copy',
           pattern: "*.{bam,out,tab}",
           saveAs: { filename -> "${meta.sample_id}/${filename}" }

    input:
        tuple val(meta), path(read1), path(read2)
        path star_index

    output:
        tuple val(meta),
              path("${meta.sample_id}.Aligned.sortedByCoord.out.bam"),
              emit: bam

        tuple val(meta),
              path("${meta.sample_id}.Log.final.out"),
              emit: log

        tuple val(meta),
              path("${meta.sample_id}.ReadsPerGene.out.tab"),
              emit: gene_counts

        tuple val(meta),
              path("${meta.sample_id}.SJ.out.tab"),
              emit: sj

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