nextflow.enable.dsl = 2

process CALL_TRANSLATED_ORFS {

    tag "${sample_id}"

    conda "conda-forge::pandas bioconda::pysam"

    publishDir "${params.outdir}/ribo/orf_calling/translated_orfs",
        mode: "copy",
        overwrite: true

    input:
    // from qc.psite_offset_qc: (sample_id, bam, bai, offsets)
    tuple val(sample_id),
          path(bam),
          path(bam_index),
          path(psite_offsets)

    // gtf passed separately (same for all samples)
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}.translated_orfs.tsv"), emit: tsv

    script:
    """
    python ${moduleDir}/call_translated_orfs.py \
        --bam ${bam} \
        --gtf ${gtf} \
        --psite_offsets ${psite_offsets} \
        --out_prefix ${sample_id}
    """
}

workflow TRANSLATED_ORFS {

    take:
    psite_offset_qc_ch
    gtf

    main:
    out = CALL_TRANSLATED_ORFS(psite_offset_qc_ch, gtf)

    emit:
    tsv = out.tsv
}