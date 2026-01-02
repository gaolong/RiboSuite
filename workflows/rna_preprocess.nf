include { FASTP } from '../modules/fastp/main.nf'

workflow RNA_PREPROCESS {

    main:

    /*
     * Parse RNA samplesheet
     * Emit: tuple(meta, R1, R2)
     */
    samples_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [ sample_id: row.sample_id ]
            tuple(
                meta,
                file(row.read1),
                file(row.read2)
            )
        }

    /*
     * fastp trimming (PE, auto-adapter, QC included)
     * Emits: tuple(meta, trimmed_R1, trimmed_R2)
     */
    FASTP(samples_ch)

    emit:
    FASTP.out.reads
    FASTP.out.html
    FASTP.out.json
}
