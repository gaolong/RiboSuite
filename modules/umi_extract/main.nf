process UMI_EXTRACT {

    tag "$sample_id"

    conda "bioconda::umi_tools=1.1.4"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.umi.fastq.gz"),
          path("${sample_id}.umi_extract.log")

    script:
    """
    set -euo pipefail

    umi_tools extract \
      --extract-method ${params.umi_extract_method} \
      --bc-pattern '${params.umi_bc_pattern}' \
      --stdin ${reads} \
      --stdout ${sample_id}.umi.fastq.gz \
      --log ${sample_id}.umi_extract.log
    """
}