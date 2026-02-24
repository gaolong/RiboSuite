process CUTADAPT_TRIM {

    tag "$sample_id"

    conda "bioconda::cutadapt=4.9 bioconda::seqtk"

    publishDir "${params.outdir}/ribo/cutadapt",
        mode: 'copy',
        saveAs: { file ->
            def fname = file instanceof Path ? file.name : file.toString()

            // always publish logs
            if (fname.endsWith('.log')) {
                return "${sample_id}/${fname}"
            }

            // publish FASTQ only if enabled
            if (fname.endsWith('.fastq.gz') && params.publish_fastq) {
                return "${sample_id}/${fname}"
            }

            // otherwise: do not publish
            return null
        }

    input:
        tuple val(sample_id),
              path(reads),
              val(adapter_5),
              val(adapter_3),
              val(bc_pattern)

    output:
        tuple val(sample_id),
              path("${sample_id}.trimmed.fastq.gz"),
              val(bc_pattern),
              emit: trimmed_reads

        path "${sample_id}.cutadapt.log", emit: cutadapt_log
        path "${sample_id}.failed.cutadapt.log", optional: true, emit: failed_log

    script:
    """
    set -euo pipefail

    FAILED_INITIAL=0

    ############################
    # Build adapter args
    ############################
    ADAPTER_ARGS=""

    if [ "${adapter_5}" != "null" ]; then
        ADAPTER_ARGS="\$ADAPTER_ARGS -g ${adapter_5}"
    fi

    if [ "${adapter_3}" != "null" ]; then
        ADAPTER_ARGS="\$ADAPTER_ARGS -a ${adapter_3}"
    fi

    ############################
    # 1. Primary trimming
    ############################
    cutadapt \\
      -j ${task.cpus} \\
      --trim-n \\
      \$ADAPTER_ARGS \\
      -a "G{10}" \\
      -m ${params.cutadapt_min_len} \\
      -M ${params.cutadapt_max_len} \\
      -o ${sample_id}.trimmed.fastq.gz \\
      ${reads} \\
      > ${sample_id}.cutadapt.log

    ############################
    # 2. Check trimming success
    ############################
    written_frac_pct=\$(
      grep -oP 'Reads written \\(passing filters\\):.*\\(\\K[0-9.]+(?=%\\))' \\
        "${sample_id}.cutadapt.log" | tail -n 1 || echo 0
    )

    adapter_frac_pct=\$(
      grep -oP 'Reads with adapters:.*\\(\\K[0-9.]+(?=%\\))' \\
        "${sample_id}.cutadapt.log" | tail -n 1 || echo 0
    )

    written_frac=\$(awk -v x="\$written_frac_pct" 'BEGIN{print x/100}')
    adapter_frac=\$(awk -v x="\$adapter_frac_pct" 'BEGIN{print x/100}')

    echo "[CUTADAPT_TRIM] ${sample_id}: written_frac=\${written_frac}, adapter_frac=\${adapter_frac}"

    if ! awk -v w="\$written_frac" -v wt="${params.cutadapt_success_frac}" \\
             -v a="\$adapter_frac" -v at="${params.cutadapt_min_adapter_frac}" \\
             'BEGIN{exit !(w >= wt && a >= at)}'; then

        echo "[CUTADAPT_TRIM] ${sample_id}: trimming failed, will attempt rescue"
        FAILED_INITIAL=1
        mv "${sample_id}.cutadapt.log" "${sample_id}.failed.cutadapt.log"
    else
        echo "[CUTADAPT_TRIM] ${sample_id}: trimming successful"
    fi

    ############################
    # 3. Fallback: adapter rescue (3′ only)
    ############################
    if [ "\$FAILED_INITIAL" -eq 1 ]; then

        seqtk sample -s42 ${reads} ${params.cutadapt_subsample_n} > sub.fastq

        best_adapter=""
        best_score=0

        adapters=(
          "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
          "AGATCGGAAGAGCACACGTCT"
          "CTGGAATTCTCGGGTGCCAAG"
          "CTGTAGGCACCATCAAT"
          "AGCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        )

        for a in "\${adapters[@]}"; do
            cutadapt \\
              -a "\$a" \\
              -O ${params.cutadapt_min_overlap} \\
              -m ${params.cutadapt_min_len} \\
              -M ${params.cutadapt_max_len} \\
              sub.fastq \\
              -o test.fastq.gz \\
              > test.log

            kept=\$(zcat test.fastq.gz | awk 'END{print NR/4}')

            if [ "\$kept" -gt "\$best_score" ]; then
                best_score="\$kept"
                best_adapter="\$a"
            fi
        done

        if [ -z "\$best_adapter" ] || [ "\$best_score" -eq 0 ]; then
            echo "[ERROR] ${sample_id}: adapter rescue failed" >&2
            exit 1
        fi

        echo "[CUTADAPT_TRIM] ${sample_id}: rescued 3' adapter = \$best_adapter"

        cutadapt \\
          -j ${task.cpus} \\
          --trim-n \\
          -a "\$best_adapter" \\
          -a "G{10}" \\
          -m ${params.cutadapt_min_len} \\
          -M ${params.cutadapt_max_len} \\
          -o ${sample_id}.trimmed.fastq.gz \\
          ${reads} \\
          >> ${sample_id}.cutadapt.log
    fi
    """
}
