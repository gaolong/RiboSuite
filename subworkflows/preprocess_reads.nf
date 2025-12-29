nextflow.enable.dsl=2

include { CUTADAPT_TRIM as CUTADAPT } from '../modules/cutadapt/main.nf'
include { UMI_EXTRACT }               from '../modules/umi_extract/main.nf'
include { FASTQC } from '../modules/fastqc/main.nf'

workflow PREPROCESS_READS {

    take:
        reads_ch   // tuple(sample_id, fastq, adapter, bc_pattern)

    main:

        /*
         * Step 0: sanity check bc_pattern
         */
        reads_checked = reads_ch.map { sid, fq, adapter, bc ->
            assert !bc || (bc.contains('(?P<') && bc.startsWith('^') && bc.endsWith('$')) :
                "Malformed bc_pattern for sample ${sid}: '${bc}'"
            tuple(sid, fq, adapter, bc)
        }

        /*
         * Step 1: adapter + polyG trimming
         */
        trimmed_ch = CUTADAPT(reads_checked).trimmed_reads
        // (sample_id, trimmed_fastq, bc_pattern)

        /*
         * Step 2: split by presence of bc_pattern
         */
        with_umi_ch = trimmed_ch.filter { sid, fq, bc -> bc }
        no_umi_ch   = trimmed_ch.filter { sid, fq, bc -> !bc }

        /*
         * Step 3: UMI extraction (SELECT FASTQ CHANNEL)
         */
        umi_fastq_ch = UMI_EXTRACT(with_umi_ch).umi_fastq
        // (sample_id, umi_fastq)

        /*
         * Step 4: pass-through for no-UMI samples
         */
        no_umi_fastq_ch = no_umi_ch.map { sid, fq, bc ->
            tuple(sid, fq)
        }

        /*
         * Step 5: merge
         */
        processed = umi_fastq_ch.mix(no_umi_fastq_ch)

        /*
        * Step 6: QC (no trimming)
        */
        qc = FASTQC(processed)

    emit:
        processed_fastq = processed
}
