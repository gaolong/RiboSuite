process FEATURECOUNTS {

    tag "${meta?.id ?: 'sample'}"
    conda "bioconda::subread=2.0.6"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf

    output:
    path "${meta?.id ?: 'sample'}.featureCounts.txt"
    path "${meta?.id ?: 'sample'}.featureCounts.txt.summary"

    script:
    def stranded = params.rna_strandedness ?: 0
    def paired   = params.rna_paired ? "-p" : ""

    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o ${meta?.id ?: 'sample'}.featureCounts.txt \\
        -t exon \\
        -g gene_id \\
        -s ${stranded} \\
        ${paired} \\
        ${bam}
    """
}
