include { FASTQC } from '../modules/fastqc/main.nf'

workflow RNA_QC_RAW {

    take:
    reads_ch   // tuple(meta, [read1, read2])

    main:
    fastqc_in_ch = reads_ch
        .flatMap { meta, reads ->
            reads.collect { fq ->
                tuple(meta, fq)
            }
        }

    FASTQC(fastqc_in_ch)

    reads_raw        = reads_ch
    fastqc_html_raw  = FASTQC.out.html
    fastqc_zip_raw   = FASTQC.out.zip

    emit:
    reads_raw
    fastqc_html_raw
    fastqc_zip_raw
}
