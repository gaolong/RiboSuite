process UMI_EXTRACT {

    tag "$sample_id"

    publishDir "${params.outdir}/preprocess/umi_extract",
               mode: 'copy',
               saveAs: { file ->
                   file instanceof Path ? "${sample_id}/${file.name}" : null
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
