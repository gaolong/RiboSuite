process MERGE_JUNCTIONS {

  tag "$sj_set_id"

  publishDir "${params.outdir}/rna/align/twopass_junctions/",
      mode: 'copy',
      pattern: "*.sjdb.chrStartEnd.tsv"

  input:
  tuple val(sj_set_id), path(sj_tabs)

  output:
  tuple val(sj_set_id), path("merged.sjdb.chrStartEnd.tsv"), emit: sjdb

  script:
  """
  awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4}' ${sj_tabs} \
    | sort -k1,1 -k2,2n -k3,3n -k4,4n \
    | uniq \
    > merged.sjdb.chrStartEnd.tsv
  """
}