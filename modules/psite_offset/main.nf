nextflow.enable.dsl = 2

process PSITE_OFFSET {

    tag "${meta.sample_id}"

    conda "bioconda::pysam=0.22.1 bioconda::samtools=1.19 conda-forge::python=3.10"

    publishDir "${params.outdir}/ribo/psite_offset",
        mode: 'copy',
        pattern: "*.tsv"

    input:
        tuple val(meta), path(bam), path(bai), path(gtf)

    output:
        tuple val(meta),
              path(bam),
              path(bai),
              path("${meta.sample_id}.psite_offsets.tsv")

    script:
    """
    python ${moduleDir}/psite_offset.py \
      --bam ${bam} \
      --gtf ${gtf} \
      --sample_id ${meta.sample_id} \
      --out ${meta.sample_id}.psite_offsets.tsv
    """
}