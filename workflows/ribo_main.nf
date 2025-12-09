nextflow.enable.dsl=2

workflow RIBO_WORKFLOW {
    Channel.empty() | view { "RiboSuite pipeline initialized." }
}