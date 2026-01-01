include { FASTQC } from '../modules/fastqc/main.nf'

workflow RNA_QC_RAW {

    take:
    reads_ch

    main:
    reads_ch
        .map { meta, reads -> reads.collect { fq -> tuple(meta, fq) } }
        .flatten()
        | FASTQC

    emit:
    reads_ch        as reads_raw
    FASTQC.out      as fastqc_raw
}
