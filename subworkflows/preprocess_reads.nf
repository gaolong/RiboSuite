nextflow.enable.dsl=2

include { UMI_EXTRACT } from '../modules/umi_extract/main.nf'
include { FASTP }      from '../modules/fastp/main.nf'

workflow PREPROCESS_READS {

    take:
        reads_ch   // tuple(sample_id, fastq)

    main:
        umi_reads = params.enable_umi \
            ? UMI_EXTRACT(reads_ch).map { sid, fq, log -> tuple(sid, fq) }
            : reads_ch

        fastp_out = FASTP(umi_reads)

    emit:
        trimmed_fastq = fastp_out[0]
}
