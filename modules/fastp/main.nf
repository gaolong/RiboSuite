nextflow.enable.dsl = 2

process FASTP {

    tag "${meta.sample_id}"

    conda "bioconda::fastp=0.23.4"

    publishDir "${params.outdir}/rna/fastp",
        mode: 'copy',
        saveAs: { file -> "${meta.sample_id}/${file}" }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.clean*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.sample_id}.fastp.json"),      emit: json
    tuple val(meta), path("${meta.sample_id}.fastp.html"),      emit: html

    script:
    def read_list = (reads instanceof List) ? reads : [reads]
    def n_reads   = read_list.size()

    if (n_reads == 2) {
        """
        fastp \
          --in1 ${read_list[0]} \
          --in2 ${read_list[1]} \
          --out1 ${meta.sample_id}.clean_1.fastq.gz \
          --out2 ${meta.sample_id}.clean_2.fastq.gz \
          --thread ${task.cpus} \
          --json ${meta.sample_id}.fastp.json \
          --html ${meta.sample_id}.fastp.html
        """
    } else {
        """
        fastp \
          --in1 ${read_list[0]} \
          --out1 ${meta.sample_id}.clean.fastq.gz \
          --thread ${task.cpus} \
          --json ${meta.sample_id}.fastp.json \
          --html ${meta.sample_id}.fastp.html
        """
    }
}