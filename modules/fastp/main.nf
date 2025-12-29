nextflow.enable.dsl=2

process FASTP {

    tag "$sample_id"

    conda "bioconda::fastp=0.23.4"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}.fastp.json", emit: fastp_json
    path "${sample_id}.fastp.html", emit: fastp_html
    path "${sample_id}.fastp.fastq.gz", emit: fastp_fastq

    script:
    """
    fastp \
      -i ${reads} \
      -o ${sample_id}.fastp.fastq.gz \
      --disable_adapter_trimming \
      --disable_quality_filtering \
      --disable_length_filtering \
      -j ${sample_id}.fastp.json \
      -h ${sample_id}.fastp.html
    """
}
