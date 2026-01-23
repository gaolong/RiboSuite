process CUTADAPT_TRIM {

    tag "$sample_id"

    conda "bioconda::cutadapt=4.9 bioconda::seqtk"

    publishDir "${params.outdir}/cutadapt",
        mode: 'copy',
        saveAs: { file -> "${sample_id}/${file}" }

    input:
        tuple val(sample_id), path(reads), val(adapter), val(bc_pattern)

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
    # 1. Primary trimming
    ############################
    cutadapt \\
      -j ${task.cpus} \\
      --trim-n \\
      -a ${adapter} \\
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
    # 3. Fallback: adapter rescue
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

        echo "[CUTADAPT_TRIM] ${sample_id}: rescued adapter = \$best_adapter"

        ############################
        # 4. Re-trim with rescued adapter
        ############################
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