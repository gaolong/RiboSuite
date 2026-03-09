nextflow.enable.dsl = 2

process BOWTIE2_FILTER {

    tag "${meta.sample_id}"

    conda "bioconda::bowtie2=2.5.2"

    /*
     * Single publishDir:
     * - logs always published
     * - FASTQ published only if params.publish_fastq
     */
    publishDir "${params.outdir}/ribo/align/bowtie2_filter",
        mode: 'copy',
        saveAs: { file ->
            def fname = file instanceof Path ? file.name : file.toString()

            // always publish logs
            if (fname.endsWith('.log')) {
                return "${meta.sample_id}/${fname}"
            }

            // publish FASTQ only if enabled
            if (fname.endsWith('.fastq.gz') && params.publish_fastq) {
                return "${meta.sample_id}/${fname}"
            }

            // otherwise: do not publish
            return null
        }

    input:
        tuple val(meta), path(reads), path(index_dir), val(index_prefix)

    output:
        tuple val(meta),
              path("${meta.sample_id}.clean.fastq.gz"),
              emit: clean_reads

        tuple val(meta),
              path("${meta.sample_id}.bowtie2.log"),
              emit: bowtie2_log

    script:
    """
    set -euo pipefail

    bowtie2 \\
      -x ${index_dir}/${index_prefix} \\
      -U ${reads} \\
      --very-sensitive \\
      -p ${task.cpus} \\
      --un ${meta.sample_id}.clean.fastq \\
      2> ${meta.sample_id}.bowtie2.log

    gzip ${meta.sample_id}.clean.fastq
    """
}