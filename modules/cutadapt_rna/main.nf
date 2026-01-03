nextflow.enable.dsl = 2

process CUTADAPT_RNA {

    tag "${sample_id}"

    conda "bioconda::cutadapt=4.9"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fq.gz"),
          path("${sample_id}_R2.trimmed.fq.gz"),
          emit: trimmed_reads

    path "${sample_id}.cutadapt.log", emit: cutadapt_log

    script:
    """
    cutadapt \
      -j ${task.cpus} \
      -o ${sample_id}_R1.trimmed.fq.gz \
      -p ${sample_id}_R2.trimmed.fq.gz \
      ${read1} ${read2} \
      > ${sample_id}.cutadapt.log
    """
}
