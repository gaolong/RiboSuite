nextflow.enable.dsl=2

process RPF_LENGTH_QC {

    tag "$sample_id"
    conda "bioconda::samtools=1.18"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.rpf_length.tsv")

    script:
    """
    # Extract read length from BAM
    samtools view $bam \
      | awk '{
          len = length(\$10);
          count[len]++
        }
        END {
          for (l in count)
            print l "\\t" count[l]
        }' \
      | sort -n > ${sample_id}.rpf_length.tsv
    """
}