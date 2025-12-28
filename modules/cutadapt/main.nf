nextflow.enable.dsl=2

process CUTADAPT_TRIM {

    tag "$sample_id"

    conda "bioconda::cutadapt=4.9"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}.cutadapt.log", emit: cutadapt_log

    script:
    """
    cutadapt \
      -a ${params.adapter_sequence} \
      -m ${params.min_length} \
      -M ${params.max_length} \
      -o ${sample_id}.trimmed.fastq.gz \
      ${reads} \
      > ${sample_id}.cutadapt.log
    """
}
