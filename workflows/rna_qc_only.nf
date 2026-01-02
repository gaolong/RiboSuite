include { RNA_QC_RAW } from '../subworkflows/rna_qc_raw.nf'

workflow RNA_QC_ONLY {

    main:
    samples_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [ sample_id: row.sample_id ]
            def reads = [ file(row.read1), file(row.read2) ]
            tuple(meta, reads)
        }

    RNA_QC_RAW(samples_ch)
}
