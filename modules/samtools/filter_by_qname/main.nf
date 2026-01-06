process SAMTOOLS_FILTER_BY_QNAME {

    tag "$sample_id"

    conda "bioconda::samtools=1.20"

    input:
        tuple val(sample_id), path(bam), path(qnames)

    output:
        tuple val(sample_id),
              path("${sample_id}.Aligned.toTranscriptome.filtered.sorted.bam"),
              path("${sample_id}.Aligned.toTranscriptome.filtered.sorted.bam.bai")

    script:
        """
        set -euo pipefail

        # 1) Filter transcriptome BAM using genome-space survivor QNAMEs
        samtools view \
          -N ${qnames} \
          -b ${bam} \
          > ${sample_id}.Aligned.toTranscriptome.filtered.bam

        # 2) Coordinate sort (REQUIRED before indexing)
        samtools sort \
          -o ${sample_id}.Aligned.toTranscriptome.filtered.sorted.bam \
          ${sample_id}.Aligned.toTranscriptome.filtered.bam

        # 3) Index sorted BAM
        samtools index \
          ${sample_id}.Aligned.toTranscriptome.filtered.sorted.bam
        """
}

