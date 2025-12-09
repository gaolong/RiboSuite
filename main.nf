nextflow.enable.dsl=2

include { RIBO_WORKFLOW } from './workflows/ribo_main.nf'

workflow {
    RIBO_WORKFLOW()
}
