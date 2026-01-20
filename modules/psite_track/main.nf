process PSITE_TRACK {

    tag "${sample_id}"

    conda "conda-forge::pandas bioconda::pysam bioconda::samtools=1.18 bioconda::bedtools=2.31.0 bioconda::ucsc-bedgraphtobigwig"

    publishDir "${params.outdir}/psite_tracks/",
        mode: 'copy',
        pattern: "*.bw"

    input:
        tuple val(sample_id), path(bam)
        path offset_tsv
        path genome_sizes

    output:
        tuple val(sample_id), path("${sample_id}.psite.len*.bw"), optional: true, emit: psite_tracks_by_len
        tuple val(sample_id), path("${sample_id}.psite.all.bw"),  emit: psite_track_all

    script:
    """
    set -euo pipefail

    # --------------------------------------------------
    # Generate P-site BEDs
    #   - per length: *.psite.lenL.bed
    #   - combined:   *.psite.all.bed
    # --------------------------------------------------
    python ${launchDir}/bin/bam_to_psite.py \
        --bam ${bam} \
        --offset ${offset_tsv} \
        --out ${sample_id}.psite \
        ${params.psite_by_length ? "--by_length" : ""}

    # --------------------------------------------------
    # Per-read-length bigWigs (auto-detect)
    # --------------------------------------------------
    if ls ${sample_id}.psite.len*.bed 1>/dev/null 2>&1; then
        for bed in ${sample_id}.psite.len*.bed; do
            len=\$(basename \$bed | sed -E 's/.*\\.len([0-9]+)\\.bed/\\1/')

            sort -k1,1 -k2,2n \$bed \
                > ${sample_id}.psite.len\${len}.sorted.bed

            bedtools genomecov \
                -i ${sample_id}.psite.len\${len}.sorted.bed \
                -g ${genome_sizes} \
                -bg > ${sample_id}.psite.len\${len}.bedGraph

            bedGraphToBigWig \
                ${sample_id}.psite.len\${len}.bedGraph \
                ${genome_sizes} \
                ${sample_id}.psite.len\${len}.bw
        done
    fi

    # --------------------------------------------------
    # Combined all-length bigWig (ALWAYS)
    # --------------------------------------------------
    sort -k1,1 -k2,2n ${sample_id}.psite.all.bed \
        > ${sample_id}.psite.all.sorted.bed

    bedtools genomecov \
        -i ${sample_id}.psite.all.sorted.bed \
        -g ${genome_sizes} \
        -bg > ${sample_id}.psite.all.bedGraph

    bedGraphToBigWig \
        ${sample_id}.psite.all.bedGraph \
        ${genome_sizes} \
        ${sample_id}.psite.all.bw
    """
}