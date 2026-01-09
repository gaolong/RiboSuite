nextflow.enable.dsl = 2

include { PSITE_TRACK } from '../modules/psite_track/main.nf'

workflow {

    Channel
        .fromPath(params.bams)
        .map { bam ->
            def sample_id = bam.baseName.replaceAll(/\.bam$/, '')
            tuple(sample_id, bam)
        }
        .set { bam_ch }

    PSITE_TRACK(
        bam_ch,
        params.psite_offset,
        params.genome_sizes
    )
}
