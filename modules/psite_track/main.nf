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
        tuple val(sample_id), path("${sample_id}.psite.all.bw"), optional: true, emit: psite_track_all
        tuple val(sample_id), path("${sample_id}.psite.pos.bw"), optional: true, emit: psite_track_pos
        tuple val(sample_id), path("${sample_id}.psite.neg.bw"), optional: true, emit: psite_track_neg
        tuple val(sample_id), path("${sample_id}.psite.pos.len*.bw"), optional: true, emit: psite_tracks_pos_by_len
        tuple val(sample_id), path("${sample_id}.psite.neg.len*.bw"), optional: true, emit: psite_tracks_neg_by_len

    script:
    // Track mode: none|all|stranded|both
    def mode = (params.psite_track_mode ?: 'both').toString()

    """
    set -euo pipefail

    # Inject Nextflow var as a literal string into bash (no unescaped \$MODE problems)
    MODE='${mode}'

    if [ "\$MODE" = "none" ]; then
      echo "PSITE_TRACK: mode=none -> skipping bigWig generation"
      exit 0
    fi

    want_all=0
    want_strand=0

    if [ "\$MODE" = "all" ] || [ "\$MODE" = "both" ]; then
      want_all=1
    fi
    if [ "\$MODE" = "stranded" ] || [ "\$MODE" = "both" ]; then
      want_strand=1
    fi

    # --------------------------------------------------
    # Generate P-site BEDs
    # --------------------------------------------------
    python ${launchDir}/bin/bam_to_psite.py \\
        --bam ${bam} \\
        --offset ${offset_tsv} \\
        --out ${sample_id}.psite \\
        ${params.psite_by_length ? "--by_length" : ""} \\
        \$([ "\$want_strand" -eq 1 ] && echo "--strand_specific" || true)

    # Helper: BED -> BW
    bed_to_bw () {
        local bed="\$1"
        local outprefix="\$2"

        # skip if bed missing/empty
        if [ ! -s "\$bed" ]; then
            echo "bed_to_bw: missing/empty bed: \$bed (skip)"
            return 0
        fi

        sort -k1,1 -k2,2n "\$bed" > "\${outprefix}.sorted.bed"

        bedtools genomecov \\
            -i "\${outprefix}.sorted.bed" \\
            -g ${genome_sizes} \\
            -bg > "\${outprefix}.bedGraph"

        bedGraphToBigWig \\
            "\${outprefix}.bedGraph" \\
            ${genome_sizes} \\
            "\${outprefix}.bw"
    }

    # --------------------------------------------------
    # Per-read-length bigWigs (combined) - only if want_all
    # --------------------------------------------------
    if [ "\$want_all" -eq 1 ]; then
      if ls ${sample_id}.psite.len*.bed 1>/dev/null 2>&1; then
          for bed in ${sample_id}.psite.len*.bed; do
              len=\$(basename "\$bed" | sed -E 's/.*\\.len([0-9]+)\\.bed/\\1/')
              bed_to_bw "\$bed" "${sample_id}.psite.len\${len}"
          done
      fi
    fi

    # --------------------------------------------------
    # Combined all-length bigWig - only if want_all
    # --------------------------------------------------
    if [ "\$want_all" -eq 1 ]; then
      bed_to_bw "${sample_id}.psite.all.bed" "${sample_id}.psite.all"
    fi

    # --------------------------------------------------
    # Strand-specific bigWigs - only if want_strand
    # --------------------------------------------------
    if [ "\$want_strand" -eq 1 ]; then
        bed_to_bw "${sample_id}.psite.pos.bed" "${sample_id}.psite.pos"
        bed_to_bw "${sample_id}.psite.neg.bed" "${sample_id}.psite.neg"

        if ls ${sample_id}.psite.pos.len*.bed 1>/dev/null 2>&1; then
            for bed in ${sample_id}.psite.pos.len*.bed; do
                len=\$(basename "\$bed" | sed -E 's/.*\\.pos\\.len([0-9]+)\\.bed/\\1/')
                bed_to_bw "\$bed" "${sample_id}.psite.pos.len\${len}"
            done
        fi

        if ls ${sample_id}.psite.neg.len*.bed 1>/dev/null 2>&1; then
            for bed in ${sample_id}.psite.neg.len*.bed; do
                len=\$(basename "\$bed" | sed -E 's/.*\\.neg\\.len([0-9]+)\\.bed/\\1/')
                bed_to_bw "\$bed" "${sample_id}.psite.neg.len\${len}"
            done
        fi
    fi
    """
}