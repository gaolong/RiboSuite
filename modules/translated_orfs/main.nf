nextflow.enable.dsl = 2

process CALL_TRANSLATED_ORFS {

    tag "${meta.sample_id}"

    conda "conda-forge::pandas bioconda::pysam"

    publishDir "${params.outdir}/ribo/orf_calling/translated_orfs",
        mode: "copy",
        overwrite: true

    input:
    // tuple(meta, bam, bai, psite_offsets, gtf)
    tuple val(meta),
          path(bam),
          path(bam_index),
          path(psite_offsets),
          path(gtf)

    output:
    tuple val(meta), path("${meta.sample_id}.translated_orfs.tsv"), emit: tsv

    script:
    """
    python ${moduleDir}/call_translated_orfs.py \
        --bam ${bam} \
        --gtf ${gtf} \
        --psite_offsets ${psite_offsets} \
        --out_prefix ${meta.sample_id}
    """
}

workflow TRANSLATED_ORFS {

    take:
    // tuple(meta, bam, bai, psite_offsets, gtf)
    psite_offset_qc_ch

    main:
    out = CALL_TRANSLATED_ORFS(psite_offset_qc_ch)

    emit:
    tsv = out.tsv
}