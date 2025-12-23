nextflow.enable.dsl=2

include { RiboSuite } from './workflows/ribosuite.nf'

workflow {

    if ( !params.reads ) {
        error "ERROR: --reads is required"
    }

    /*
     * Convert params.reads into a channel of:
     * tuple(sample_id, fastq)
     */
    reads_ch = Channel.fromPath(params.reads)
        .map { file ->
            tuple(
                file.baseName.replaceAll(/\.fastq(\.gz)?$/, ''),
                file
            )
        }

    RiboSuite(reads_ch)
}