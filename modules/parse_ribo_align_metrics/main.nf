nextflow.enable.dsl = 2

process PARSE_RIBO_ALIGN_METRICS {

    tag { meta.sample_id }

    input:
        tuple val(meta), path(log_final)

    output:
        path "${meta.sample_id}.ribo_align_metrics.tsv", emit: metrics

    script:
    """
    set -euo pipefail

    unique_reads=\$(grep -F "Uniquely mapped reads number" "${log_final}" | tail -n 1 | cut -d'|' -f2 | tr -d '[:space:]')
    unique_pct=\$(grep -F "Uniquely mapped reads %" "${log_final}" | tail -n 1 | cut -d'|' -f2 | tr -d '[:space:]%' )

    [ -n "\$unique_reads" ] || unique_reads="NA"
    [ -n "\$unique_pct" ]   || unique_pct="NA"

    cat > ${meta.sample_id}.ribo_align_metrics.tsv <<EOF
sample_id\tunique_mapped_reads\tunique_mapped_pct
${meta.sample_id}\t\$unique_reads\t\$unique_pct
EOF
    """
}