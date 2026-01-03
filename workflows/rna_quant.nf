include { RNA_PREPROCESS }  from './rna_preprocess.nf'
include { ALIGN_RNA }       from '../subworkflows/align_rna.nf'
include { FEATURECOUNTS }  from '../modules/featurecounts/main.nf'

workflow RNA_QUANT {

    main:

    /*
     * Step 1: RNA preprocessing
     */
    RNA_PREPROCESS()

    trimmed_reads_ch = RNA_PREPROCESS.out[0]

    /*
     * Step 2: STAR alignment
     */
    ALIGN_RNA(trimmed_reads_ch)

    /*
     * ALIGN_RNA emits: (meta, bam, bai)
     */
    aligned_bam_ch = ALIGN_RNA.out[0]

    /*
     * Step 3: featureCounts
     */
    FEATURECOUNTS(
        aligned_bam_ch,
        params.gtf
    )
}
