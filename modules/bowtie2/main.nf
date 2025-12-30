nextflow.enable.dsl=2

process BOWTIE2_FILTER {

    tag "$sample_id"
    conda "bioconda::bowtie2=2.5.2"

    input:
        tuple val(sample_id), path(reads)
        path index_dir
        val  index_prefix

    output:
        tuple val(sample_id),
              path("${sample_id}.clean.fastq.gz"), emit: clean_reads
        path "${sample_id}.bowtie2.log", emit: bowtie2_log

    script:
        """
        bowtie2 \
          -x ${index_dir}/${index_prefix} \
          -U ${reads} \
          --very-sensitive \
          -p ${task.cpus} \
          --un ${sample_id}.clean.fastq \
          2> ${sample_id}.bowtie2.log

        gzip ${sample_id}.clean.fastq
        """
}
