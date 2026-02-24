process FASTQC {

    tag "$sample_id"

    conda "bioconda::fastqc=0.12.1"

    publishDir "${params.outdir}/ribo/qc/fastqc",
        mode: 'copy',
        pattern: "*_fastqc.{html,zip}",
        saveAs: { file -> "${sample_id}/${file}" }

    input:
        tuple val(sample_id), path(reads)

    output:
        path "*_fastqc.html", emit: html
        path "*_fastqc.zip",  emit: zip

    script:
    """
    set -euo pipefail
    fastqc ${reads} --outdir .
    """
}
