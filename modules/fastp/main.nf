nextflow.enable.dsl = 2

process FASTP {

    tag "${meta.sample_id}"

    conda "bioconda::fastp=0.23.4"

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta),
          path("${meta.sample_id}_R1.fastp.fq.gz"),
          path("${meta.sample_id}_R2.fastp.fq.gz"),
          emit: reads

    path "${meta.sample_id}.fastp.json", emit: json
    path "${meta.sample_id}.fastp.html", emit: html

    script:
    """
    fastp \
      --in1 ${read1} \
      --in2 ${read2} \
      --out1 ${meta.sample_id}_R1.fastp.fq.gz \
      --out2 ${meta.sample_id}_R2.fastp.fq.gz \
      --detect_adapter_for_pe \
      --thread ${task.cpus} \
      --json ${meta.sample_id}.fastp.json \
      --html ${meta.sample_id}.fastp.html
    """
}
