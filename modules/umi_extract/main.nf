process UMI_EXTRACT {

    tag "$sample_id"

    conda "bioconda::umi_tools=1.1.4"

    /*
     * Single publishDir:
     * - logs always published
     * - FASTQ published only if params.publish_fastq
     */
    publishDir "${params.outdir}/ribo/preprocess/umi_extract",
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
        tuple val(sample_id), path(reads), val(bc_pattern)

    output:
        tuple val(sample_id),
              path("${sample_id}.umi.fastq.gz"),
              emit: umi_fastq

        path "${sample_id}.umi_extract.log", emit: umi_log

    script:
    """
    set -euo pipefail

    umi_tools extract \
      --extract-method=regex \
      --bc-pattern='${bc_pattern}' \
      --stdin ${reads} \
      --stdout ${sample_id}.umi.fastq.gz \
      --log ${sample_id}.umi_extract.log
    """
}
