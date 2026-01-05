nextflow.enable.dsl=2

process CUTADAPT_TRIM {

    tag "$sample_id"

    conda "bioconda::cutadapt=4.9"

    input:
    tuple val(sample_id), path(reads), val(adapter), val(bc_pattern)

    output:
    tuple val(sample_id),
          path("${sample_id}.trimmed.fastq.gz"),
          val(bc_pattern),
          emit: trimmed_reads

    path "${sample_id}.cutadapt.log", emit: cutadapt_log

    script:
    """
    cutadapt \
      -j ${task.cpus} \
      -a ${adapter} \
      -a "G{10}" \
      -o ${sample_id}.trimmed.fastq.gz \
      ${reads} \
      > ${sample_id}.cutadapt.log
    """
}
