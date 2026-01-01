include { RNA_QC_RAW } from '../subworkflows/rna_qc_raw.nf'

workflow RNA_QC_ONLY {

    main:
    samples_ch = Channel.fromPath(params.samplesheet)
        | parseSamplesheet   // reuse your existing logic

    RNA_QC_RAW(samples_ch)
}
