process FASTP {
    tag "$sample_id"
    conda "bioconda::fastp=0.23.4"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), 
          path("${sample_id}.trimmed.fastq.gz"), 
          path("${sample_id}.fastp.json"), 
          path("${sample_id}.fastp.html")

    script:
    """
    fastp \
        -i $reads \
        -o ${sample_id}.trimmed.fastq.gz \
        --json ${sample_id}.fastp.json \
        --html ${sample_id}.fastp.html
    """
}