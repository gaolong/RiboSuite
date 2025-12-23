nextflow.enable.dsl=2

process FASTP {

    tag "$sample_id"

    conda "bioconda::fastp=0.23.4"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}.fastp.json", emit: fastp_json
    path "${sample_id}.fastp.html", emit: fastp_html

    script:
    """
    fastp \
      -i ${reads} \
      -o ${sample_id}.trimmed.fastq.gz \
      -j ${sample_id}.fastp.json \
      -h ${sample_id}.fastp.html
    """
}