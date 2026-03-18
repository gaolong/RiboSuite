nextflow.enable.dsl = 2

process RNA_TRACK {

  tag "${meta.sample_id}"

  conda "bioconda::samtools=1.18 bioconda::bedtools=2.31.0 bioconda::ucsc-bedgraphtobigwig"

  publishDir "${params.outdir}/rna/tracks/",
      mode: 'copy',
      pattern: "*.bw"

  input:
    tuple val(meta), path(bam), path(genome_sizes)

  output:
    tuple val(meta), path("${meta.sample_id}.rna.all.bw"), optional: true, emit: rna_track_all
    tuple val(meta), path("${meta.sample_id}.rna.pos.bw"), optional: true, emit: rna_track_pos
    tuple val(meta), path("${meta.sample_id}.rna.neg.bw"), optional: true, emit: rna_track_neg

  script:
    def sample_id   = meta.sample_id
    def mode        = (params.rna_track_mode ?: 'both').toString()
    def stranded    = (params.rna_track_strandedness != null) ? params.rna_track_strandedness as Integer : 0
    def is_paired   = meta?.containsKey('is_paired') ? meta.is_paired : false
    def paired_flag = is_paired ? "-pc" : ""

    """
    set -euo pipefail

    MODE="${mode}"
    STRANDED="${stranded}"

    want_all=0
    want_strand=0

    if [ "\$MODE" = "all" ] || [ "\$MODE" = "both" ]; then
      want_all=1
    fi
    if [ "\$MODE" = "stranded" ] || [ "\$MODE" = "both" ]; then
      want_strand=1
    fi

    # If library is unstranded (0), disable stranded output
    if [ "\$STRANDED" -eq 0 ]; then
      want_strand=0
    fi

    # --------------------------------------------------
    # All-strand bigWig
    # --------------------------------------------------
    if [ "\$want_all" -eq 1 ]; then
      bedtools genomecov ${paired_flag} -ibam ${bam} -bg -split \\
        | sort -k1,1 -k2,2n \\
        > ${sample_id}.rna.all.bedGraph

      bedGraphToBigWig \\
        ${sample_id}.rna.all.bedGraph \\
        ${genome_sizes} \\
        ${sample_id}.rna.all.bw
    fi

    # --------------------------------------------------
    # Strand-specific bigWigs
    # --------------------------------------------------
    if [ "\$want_strand" -eq 1 ]; then

      # featureCounts convention:
      # 1 = same strand
      # 2 = reverse strand
      if [ "\$STRANDED" -eq 1 ]; then
        POS_STRAND="+"
        NEG_STRAND="-"
      else
        POS_STRAND="-"
        NEG_STRAND="+"
      fi

      # positive
      bedtools genomecov ${paired_flag} -ibam ${bam} -bg -split -strand "\${POS_STRAND}" \\
        | sort -k1,1 -k2,2n \\
        > ${sample_id}.rna.pos.bedGraph

      bedGraphToBigWig \\
        ${sample_id}.rna.pos.bedGraph \\
        ${genome_sizes} \\
        ${sample_id}.rna.pos.bw

      # negative
      bedtools genomecov ${paired_flag} -ibam ${bam} -bg -split -strand "\${NEG_STRAND}" \\
        | sort -k1,1 -k2,2n \\
        > ${sample_id}.rna.neg.bedGraph

      bedGraphToBigWig \\
        ${sample_id}.rna.neg.bedGraph \\
        ${genome_sizes} \\
        ${sample_id}.rna.neg.bw
    fi
    """
}