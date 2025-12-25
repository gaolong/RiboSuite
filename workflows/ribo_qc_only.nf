nextflow.enable.dsl=2

include { RIBO_QC_BASIC } from '../subworkflows/ribo_qc_basic.nf'

workflow {

    bam_file = file(params.bam)
    bai_file = file("${params.bam}.bai")
    gtf_file = file(params.gtf)

    if( !bai_file.exists() ) {
        error "BAM index not found: ${bai_file}"
    }

    aligned_ch = Channel.of(
        tuple(
            bam_file.baseName.replace('.Aligned.sortedByCoord.out',''),
            bam_file,
            bai_file
        )
    )

    RIBO_QC_BASIC(aligned_ch, gtf_file)
}
