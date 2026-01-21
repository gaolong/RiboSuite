process CUTADAPT_TRIM {

    tag "$sample_id"

    conda "bioconda::cutadapt=4.9 bioconda::seqtk"

    publishDir "${params.outdir}/cutadapt", mode: 'copy'

    input:
        tuple val(sample_id), path(reads), val(adapter), val(bc_pattern)

    output:
        tuple val(sample_id),
              path("${sample_id}.trimmed.fastq.gz"),
              val(bc_pattern),
              emit: trimmed_reads

        path "${sample_id}.cutadapt.log", emit: cutadapt_log

    script:
    """
    set -euo pipefail

    ############################
    # 1. Primary trimming
    ############################
    cutadapt \\
      -j ${task.cpus} \\
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
    written_frac=\$(grep -oP "Reads written \\(passing filters\\):.*\\(\\K[0-9.]+(?=%\\))" \
                  ${sample_id}.cutadapt.log | tail -n 1)

    echo "[CUTADAPT_TRIM] ${sample_id}: written fraction = \$written_frac%"

    if awk -v f=\$written_frac -v th=${params.cutadapt_success_frac} \
          'BEGIN{exit !(f >= th)}'; then
        echo "[CUTADAPT_TRIM] ${sample_id}: trimming successful, proceeding"
        exit 0
    else
        echo "[CUTADAPT_TRIM] ${sample_id}: trimming failed, renaming initial log"
        mv ${sample_id}.cutadapt.log ${sample_id}.cutadapt.failed.log
    fi


    ############################
    # 3. Fallback: adapter rescue
    ############################
    echo "[CUTADAPT_TRIM] ${sample_id}: trimming failed, running adapter rescue"

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
      -a "\$best_adapter" \\
      -a "G{10}" \\
      -m ${params.cutadapt_min_len} \\
      -M ${params.cutadapt_max_len} \\
      -o ${sample_id}.trimmed.fastq.gz \\
      ${reads} \\
      >> ${sample_id}.cutadapt.log
    """
}
