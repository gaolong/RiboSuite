process TRANSLATED_ORFS {

    tag "${sample_id}"

    publishDir "${params.outdir}/orf_calling/translated_orfs",
        mode: "copy"

    input:
    tuple val(sample_id),
          path(bam),
          path(bam_index),
          path(gtf),
          path(psite_offsets)

    output:
    path("${sample_id}.translated_orfs.tsv")

    script:
    """
    call_translated_orfs.py \
        --bam ${bam} \
        --gtf ${gtf} \
        --psite_offsets ${psite_offsets} \
        --out_prefix ${sample_id}
    """
}
