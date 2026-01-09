process PSITE_TRACK {

    tag "${sample_id}"

    conda "conda-forge::pandas bioconda::pysam bioconda::samtools=1.18 bioconda::bedtools=2.31.0 bioconda::ucsc-bedgraphtobigwig"

    input:
        tuple val(sample_id), path(bam)
        path offset_tsv
        path genome_sizes

    output:
        tuple val(sample_id), path("${sample_id}.psite.len*.bw"), emit: psite_tracks_by_len

    script:
    """
    set -euo pipefail

    # --------------------------------------------------
    # Generate per-read-length P-site BEDs
    # --------------------------------------------------
    python ${launchDir}/bin/bam_to_psite.py \
        --bam ${bam} \
        --offset ${offset_tsv} \
        --out ${sample_id}.psite

    # --------------------------------------------------
    # Convert each BED to bigWig
    # --------------------------------------------------
    for bed in ${sample_id}.psite.len*.bed; do
        len=\$(basename \$bed | sed -E 's/.*\\.len([0-9]+)\\.bed/\\1/')

        sort -k1,1 -k2,2n \$bed > ${sample_id}.psite.len\${len}.sorted.bed

        bedtools genomecov \
            -i ${sample_id}.psite.len\${len}.sorted.bed \
            -g ${genome_sizes} \
            -bg > ${sample_id}.psite.len\${len}.bedGraph

        bedGraphToBigWig \
            ${sample_id}.psite.len\${len}.bedGraph \
            ${genome_sizes} \
            ${sample_id}.psite.len\${len}.bw
    done
    """
}
