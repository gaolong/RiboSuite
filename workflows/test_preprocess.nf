nextflow.enable.dsl=2

include { PREPROCESS_READS } from '../subworkflows/preprocess_reads.nf'

workflow TEST_PREPROCESS {

    main:
        Channel
            .fromPath(params.samples)
            .splitCsv(header:true, sep:'\t')
            .map { row ->
                tuple(
                    row.sample_id,
                    file(row.fastq),
                    row.adapter,
                    row.bc_pattern
                )
            }
            .set { reads_ch }

        preprocessed = PREPROCESS_READS(reads_ch)

    emit:
        preprocessed
}
