process PSITE_TRACK {

    tag "${sample_id}"

    conda "conda-forge::pandas bioconda::pysam bioconda::samtools=1.18 bioconda::bedtools=2.31.0 bioconda::ucsc-bedgraphtobigwig"

    publishDir "${params.outdir}/ribo/psite_tracks/",
        mode: 'copy',
        pattern: "*.bw"

    input:
        tuple val(sample_id), path(bam)
        path offset_tsv
        path genome_sizes

    output:
        tuple val(sample_id), path("${sample_id}.psite.len*.bw"), optional: true, emit: psite_tracks_by_len
        tuple val(sample_id), path("${sample_id}.psite.all.bw"),  emit: psite_track_all

        // NEW (optional)
        tuple val(sample_id), path("${sample_id}.psite.pos.bw"), optional: true, emit: psite_track_pos
        tuple val(sample_id), path("${sample_id}.psite.neg.bw"), optional: true, emit: psite_track_neg
        tuple val(sample_id), path("${sample_id}.psite.pos.len*.bw"), optional: true, emit: psite_tracks_pos_by_len
        tuple val(sample_id), path("${sample_id}.psite.neg.len*.bw"), optional: true, emit: psite_tracks_neg_by_len

    script:
    """
    set -euo pipefail

    # --------------------------------------------------
    # Generate P-site BEDs
    #   - per length: *.psite.lenL.bed
    #   - combined:   *.psite.all.bed
    #   - optional strand-specific:
    #       *.psite.pos.bed / *.psite.neg.bed
    #       *.psite.pos.lenL.bed / *.psite.neg.lenL.bed
    # --------------------------------------------------
    python ${launchDir}/bin/bam_to_psite.py \
        --bam ${bam} \
        --offset ${offset_tsv} \
        --out ${sample_id}.psite \
        ${params.psite_by_length ? "--by_length" : ""} \
        ${params.psite_strand_specific ? "--strand_specific" : ""}

    # Helper: BED -> BW
    bed_to_bw () {
        local bed=\$1
        local outprefix=\$2

        sort -k1,1 -k2,2n \$bed > \${outprefix}.sorted.bed

        bedtools genomecov \
            -i \${outprefix}.sorted.bed \
            -g ${genome_sizes} \
            -bg > \${outprefix}.bedGraph

        bedGraphToBigWig \
            \${outprefix}.bedGraph \
            ${genome_sizes} \
            \${outprefix}.bw
    }

    # --------------------------------------------------
    # Per-read-length bigWigs (combined)
    # --------------------------------------------------
    if ls ${sample_id}.psite.len*.bed 1>/dev/null 2>&1; then
        for bed in ${sample_id}.psite.len*.bed; do
            len=\$(basename \$bed | sed -E 's/.*\\.len([0-9]+)\\.bed/\\1/')
            bed_to_bw "\$bed" "${sample_id}.psite.len\${len}"
        done
    fi

    # --------------------------------------------------
    # Combined all-length bigWig (ALWAYS)
    # --------------------------------------------------
    bed_to_bw "${sample_id}.psite.all.bed" "${sample_id}.psite.all"

    # --------------------------------------------------
    # Strand-specific bigWigs (optional)
    # --------------------------------------------------
    if ${params.psite_strand_specific ? "true" : "false"}; then
        # all-length pos/neg
        if [ -s "${sample_id}.psite.pos.bed" ]; then
            bed_to_bw "${sample_id}.psite.pos.bed" "${sample_id}.psite.pos"
        fi
        if [ -s "${sample_id}.psite.neg.bed" ]; then
            bed_to_bw "${sample_id}.psite.neg.bed" "${sample_id}.psite.neg"
        fi

        # per-length pos/neg
        if ls ${sample_id}.psite.pos.len*.bed 1>/dev/null 2>&1; then
            for bed in ${sample_id}.psite.pos.len*.bed; do
                len=\$(basename \$bed | sed -E 's/.*\\.pos\\.len([0-9]+)\\.bed/\\1/')
                bed_to_bw "\$bed" "${sample_id}.psite.pos.len\${len}"
            done
        fi

        if ls ${sample_id}.psite.neg.len*.bed 1>/dev/null 2>&1; then
            for bed in ${sample_id}.psite.neg.len*.bed; do
                len=\$(basename \$bed | sed -E 's/.*\\.neg\\.len([0-9]+)\\.bed/\\1/')
                bed_to_bw "\$bed" "${sample_id}.psite.neg.len\${len}"
            done
        fi
    fi
    """
}