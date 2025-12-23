nextflow.enable.dsl=2

include { RiboSuite } from './workflows/ribosuite.nf'

workflow {

    Channel
        .fromPath(params.reads)
        .map { fastq ->
            def sample_id = fastq.baseName.replaceFirst(/\.fastq(\.gz)?$/, '')
            tuple(sample_id, fastq)
        }
        | RiboSuite
}