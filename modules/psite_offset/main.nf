process PSITE_OFFSET {

    tag "$sample_id"
    conda "bioconda::pysam=0.22.1 bioconda::samtools=1.19 conda-forge::python=3.10"

    input:
        tuple val(sample_id), path(bam), path(bai)
        path gtf

    output:
        tuple val(sample_id),
              path(bam),
              path(bai),
              path("${sample_id}.psite_offsets.tsv")

    script:
    """
    python ${moduleDir}/psite_offset.py \
      --bam ${bam} \
      --gtf ${gtf} \
      --sample_id ${sample_id} \
      --out ${sample_id}.psite_offsets.tsv
    """
}
