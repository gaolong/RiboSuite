nextflow.enable.dsl=2

include { RIBO_QC_BASIC } from '../subworkflows/ribo_qc_basic.nf'

workflow {

    aligned_ch = Channel.of(
        tuple(
            "SRR1173905",
            file(params.bam)
        )
    )

    RIBO_QC_BASIC(aligned_ch, file(params.gtf))
}
