include { RNA_PREPROCESS } from './rna_preprocess.nf'
include { ALIGN_RNA }      from '../subworkflows/align_rna.nf'

workflow RNA_QUANT {

    main:

    /*
     * Run RNA preprocessing (fastp)
     * RNA_PREPROCESS emits:
     *   - reads
     *   - html
     *   - json
     */
    RNA_PREPROCESS()

    /*
     * Extract ONLY the trimmed reads channel
     */
    trimmed_reads_ch = RNA_PREPROCESS.out[0]
    // or equivalently:
    // trimmed_reads_ch = RNA_PREPROCESS.out.reads

    /*
     * Align RNA-seq reads with STAR
     */
    ALIGN_RNA(trimmed_reads_ch)
}
