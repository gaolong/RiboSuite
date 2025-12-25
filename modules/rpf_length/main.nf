process RPF_LENGTH_QC {

    tag "$sample_id"
    conda "bioconda::samtools=1.19"

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.rpf_length.tsv")

    shell:
    '''
    samtools view !{bam} \
      | awk '{ print length($10) }' \
      | sort | uniq -c | sort -n \
      > !{sample_id}.rpf_length.tsv
    '''
}