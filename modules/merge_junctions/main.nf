nextflow.enable.dsl = 2

process MERGE_JUNCTIONS {

  tag "$sj_set_id"

  input:
  tuple val(sj_set_id), path(sj_tabs)

  output:
  tuple val(sj_set_id), path("${sj_set_id}.sjdb.chrStartEnd.tsv"), emit: sjdb

  script:
  """
  # SJ.out.tab: col1 chr, col2 intronStart, col3 intronEnd, col4 strand
  cat ${sj_tabs} \
    | awk 'BEGIN{OFS="\\t"} {print $1,$2,$3,$4}' \
    | sort -k1,1 -k2,2n -k3,3n -k4,4n \
    | uniq \
    > ${sj_set_id}.sjdb.chrStartEnd.tsv
  """
}