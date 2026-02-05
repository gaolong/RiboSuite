process UMI_EXTRACT {

    tag "$sample_id"

    // Always publish logs
    publishDir "${params.outdir}/preprocess/umi_extract/logs",
        mode: 'copy',
        saveAs: { file ->
            file.name.endsWith('.log')
                ? "${sample_id}/${file.name}"
                : null
        }

    // Publish FASTQ only if explicitly requested
    publishDir {
        if (params.publish_fastq) {
            path "${params.outdir}/preprocess/umi_extract/fastq"
            mode 'copy'
            saveAs: { file ->
                file.name.endsWith('.fastq.gz')
                    ? "${sample_id}/${file.name}"
                    : null
            }
        }
    }

    conda "bioconda::umi_tools=1.1.4"

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
