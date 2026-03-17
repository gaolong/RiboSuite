nextflow.enable.dsl = 2

include { CUTADAPT_TRIM as CUTADAPT } from '../modules/cutadapt/main.nf'
include { UMI_EXTRACT }               from '../modules/umi_extract/main.nf'
include { FASTQC }                    from '../modules/fastqc/main.nf'
include { BOWTIE2_FILTER }            from '../modules/bowtie2/main.nf'

workflow PREPROCESS_READS {

    take:
        // tuple(meta, fastq, adapter_5, adapter_3, bc_pattern, has_umi, contam_index)
        reads_ch

    main:
        /*
         * 0. Sanity check bc_pattern (only if has_umi)
         */
        checked = reads_ch.map { meta, fq, adapter_5, adapter_3, bc, has_umi, contam_index ->
            if (has_umi) {
                assert bc && bc.contains('(?P<') && bc.startsWith('^') && bc.endsWith('$') :
                    "Malformed bc_pattern for sample ${meta.sample_id}: '${bc}'"
            }
            tuple(meta, fq, adapter_5, adapter_3, bc, has_umi, contam_index)
        }

        /*
         * 1. Adapter trimming
         * CUTADAPT expects: (sample_id, fq, adapter_5, adapter_3, bc_pattern)
         */
        cutadapt_in = checked.map { meta, fq, adapter_5, adapter_3, bc, has_umi, contam_index ->
            tuple(meta.sample_id.toString(), fq, adapter_5, adapter_3, bc)
        }

        trimmed0 = CUTADAPT(cutadapt_in).trimmed_reads
        // (sid, trimmed_fastq, bc_pattern)

        /*
         * 2. Re-attach meta, has_umi, contam_index by sample_id
         * Result: (sid, trimmed_fastq, bc_pattern, meta, has_umi, contam_index)
         */
        meta_info = checked.map { meta, fq, adapter_5, adapter_3, bc, has_umi, contam_index ->
            tuple(meta.sample_id.toString(), meta, has_umi, contam_index)
        }

        trimmed = trimmed0.join(meta_info, by: 0)

        /*
         * 3. Split by UMI presence
         */
        with_umi    = trimmed.filter  { sid, fq, bc, meta, has_umi, contam_index -> has_umi }
        without_umi = trimmed.filter  { sid, fq, bc, meta, has_umi, contam_index -> !has_umi }

        /*
         * 4a. UMI extraction
         * UMI_EXTRACT expects: (sample_id, fq, bc_pattern)
         */
        umi_res = UMI_EXTRACT(
            with_umi.map { sid, fq, bc, meta, has_umi, contam_index ->
                tuple(sid, fq, bc)
            }
        )
        umi_fastq_ch = umi_res.umi_fastq   // (sid, umi.fastq.gz)

        /*
         * 4b. Non-UMI passthrough
         */
        no_umi_fastq_ch = without_umi.map { sid, fq, bc, meta, has_umi, contam_index ->
            tuple(sid, fq)
        }

        /*
         * 5. Pre-bowtie processed0
         * Result: (meta, fq, has_umi, contam_index)
         */
        umi_processed = umi_fastq_ch
            .join(
                with_umi.map { sid, fq, bc, meta, has_umi, contam_index ->
                    tuple(sid, meta, contam_index)
                },
                by: 0
            )
            .map { sid, fq, meta, contam_index ->
                tuple(meta, fq, true, contam_index)
            }

        no_umi_processed = no_umi_fastq_ch
            .join(
                without_umi.map { sid, fq, bc, meta, has_umi, contam_index ->
                    tuple(sid, meta, contam_index)
                },
                by: 0
            )
            .map { sid, fq, meta, contam_index ->
                tuple(meta, fq, false, contam_index)
            }

        processed0 = umi_processed.mix(no_umi_processed)

        /*
         * 6. Contaminant filtering
         * BOWTIE2_FILTER input:
         *   tuple(meta, fq, index_dir, index_prefix)
         */
        bowtie_in = processed0.map { meta, fq, has_umi, contam_index ->
            def contam_file = file(contam_index)
            tuple(
                meta,
                fq,
                file(contam_file.parent),
                contam_file.getName()
            )
        }

        bowtie_clean = BOWTIE2_FILTER(bowtie_in).clean_reads
        // (meta, clean_fastq)

        /*
         * 7. Re-attach has_umi after bowtie
         * Result: (meta, clean_fastq, has_umi)
         */
        processed = bowtie_clean
            .join(
                processed0.map { meta, fq, has_umi, contam_index ->
                    tuple(meta, has_umi)
                },
                by: 0
            )
            .map { meta, fq, has_umi ->
                tuple(meta, fq, has_umi)
            }

        /*
         * 8. QC (side-effect only)
         */
        FASTQC(
            processed.map { meta, fq, has_umi ->
                tuple(meta.sample_id.toString(), fq)
            }
        )

    emit:
        processed
}