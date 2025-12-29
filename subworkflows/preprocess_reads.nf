nextflow.enable.dsl=2

include { CUTADAPT_TRIM as CUTADAPT } from '../modules/cutadapt/main.nf'
include { UMI_EXTRACT }               from '../modules/umi_extract/main.nf'
include { FASTQC }                    from '../modules/fastqc/main.nf'

workflow PREPROCESS_READS {

    take:
        reads_ch   // (sample_id, fastq, adapter, bc_pattern, has_umi)

    main:
        /*
         * 0. Sanity check bc_pattern (only if has_umi)
         */
        checked = reads_ch.map { sid, fq, adapter, bc, has_umi ->
            if (has_umi) {
                assert bc && bc.contains('(?P<') && bc.startsWith('^') && bc.endsWith('$') :
                    "Malformed bc_pattern for sample ${sid}: '${bc}'"
            }
            tuple(sid, fq, adapter, bc, has_umi)
        }

        /*
         * 1. Adapter trimming
         * IMPORTANT: feed CUTADAPT only what it declares (likely 4 fields)
         */
        cutadapt_in = checked.map { sid, fq, adapter, bc, has_umi ->
            tuple(sid, fq, adapter, bc)
        }

        // CUTADAPT output assumed: (sid, trimmed_fastq, bc_pattern)
        trimmed0 = CUTADAPT(cutadapt_in).trimmed_reads

        /*
         * 2. Re-attach has_umi flag by join on sample_id
         */
        meta_has_umi = checked.map { sid, fq, adapter, bc, has_umi ->
            tuple(sid, has_umi)
        }

        // join gives: (sid, trimmed_fastq, bc_pattern, has_umi)
        trimmed = trimmed0.join(meta_has_umi)

        /*
         * 3. Split by UMI presence
         */
        with_umi    = trimmed.filter  { sid, fq, bc, has_umi -> has_umi }
        without_umi = trimmed.filter  { sid, fq, bc, has_umi -> !has_umi }

        /*
         * 4a. UMI extraction (only if has_umi)
         *
         * UMI_EXTRACT returns a *module output object* (multi-channel),
         * so grab the channel explicitly via .umi_fastq
         */
        umi_res = UMI_EXTRACT(
            with_umi.map { sid, fq, bc, has_umi -> tuple(sid, fq, bc) }
        )
        umi_fastq_ch = umi_res.umi_fastq   // (sid, umi.fastq.gz)

        /*
         * 4b. Pass-through for non-UMI samples
         */
        no_umi_fastq_ch = without_umi.map { sid, fq, bc, has_umi -> tuple(sid, fq) }

        /*
         * 5. Canonical output: (sample_id, clean_fastq, has_umi)
         */
        processed = umi_fastq_ch
            .map { sid, fq -> tuple(sid, fq, true) }
            .mix(
                no_umi_fastq_ch.map { sid, fq -> tuple(sid, fq, false) }
            )

        /*
         * 6. QC (side-effect only)
         */
        FASTQC(processed.map { sid, fq, has_umi -> tuple(sid, fq) })

    emit:
        processed
}
