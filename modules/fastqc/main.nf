process FASTQC {

    tag "$sample_id"

    conda "bioconda::fastqc=0.12.1"

    publishDir "${params.outdir}/qc/fastqc",
           mode: 'copy',
           saveAs: { file ->
               file instanceof Path ? "${sample_id}/${file.name}" : null
           }

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip

    script:
    """
    fastqc ${reads} --outdir .
    """
}