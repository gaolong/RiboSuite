nextflow.enable.dsl = 2

process BOWTIE2_FILTER {

    tag "$sample_id"

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
                return "${sample_id}/${fname}"
            }

            // publish FASTQ only if enabled
            if (fname.endsWith('.fastq.gz') && params.publish_fastq) {
                return "${sample_id}/${fname}"
            }

            // otherwise: do not publish
            return null
        }

    input:
        tuple val(sample_id), path(reads)
        path index_dir
        val  index_prefix

    output:
        tuple val(sample_id),
              path("${sample_id}.clean.fastq.gz"),
              emit: clean_reads

        path "${sample_id}.bowtie2.log", emit: bowtie2_log

    script:
    """
    set -euo pipefail

    bowtie2 \\
      -x ${index_dir}/${index_prefix} \\
      -U ${reads} \\
      --very-sensitive \\
      -p ${task.cpus} \\
      --un ${sample_id}.clean.fastq \\
      2> ${sample_id}.bowtie2.log

    gzip ${sample_id}.clean.fastq
    """
}
