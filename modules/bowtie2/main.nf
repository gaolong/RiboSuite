process BOWTIE2_FILTER {

    tag "$sample_id"
    cpus 2
    conda "bioconda::bowtie2=2.5.2"

    input:
    tuple val(sample_id), path(reads)
    val index_prefix

    output:
    tuple val(sample_id), path("${sample_id}.clean.fastq.gz"), emit: clean_reads
    path "${sample_id}.bowtie2.log", emit: bowtie2_log

    script:
    """
    bowtie2 \
        -x ${index_prefix} \
        -U ${reads} \
        --very-sensitive \
        -p ${task.cpus} \
        --un-gz ${sample_id}.clean.fastq.gz \
        2> ${sample_id}.bowtie2.log
    """
}