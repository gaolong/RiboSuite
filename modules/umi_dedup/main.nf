process UMI_DEDUP {

    tag "$sample_id"

    conda "bioconda::umi_tools=1.1.4 bioconda::samtools=1.20"

    publishDir "${params.outdir}/ribo/align/umi_dedup",
        mode: 'copy',
        saveAs: { file ->
            file.toString().endsWith('.umi_survivors.qnames.txt')
                ? null
                : "${sample_id}/${file}"
        }

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        tuple val(sample_id),
              path("${sample_id}.umi_dedup.bam"),
              path("${sample_id}.umi_dedup.bam.bai"),
              path("${sample_id}.umi_dedup.log"),
              path("${sample_id}.umi_survivors.qnames.txt")

    script:
    """
    set -euo pipefail

    # 1) UMI deduplication in removal mode (genome space)
    umi_tools dedup \
        --stdin ${bam} \
        --stdout ${sample_id}.umi_dedup.bam \
        --log ${sample_id}.umi_dedup.log

    # 2) Index deduplicated BAM
    samtools index ${sample_id}.umi_dedup.bam

    # 3) Survivor QNAMEs (used internally only)
    samtools view ${sample_id}.umi_dedup.bam \
        | cut -f1 \
        | sort -u \
        > ${sample_id}.umi_survivors.qnames.txt
    """
}
