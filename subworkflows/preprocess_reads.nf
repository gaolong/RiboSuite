nextflow.enable.dsl = 2

include { CUTADAPT_TRIM as CUTADAPT } from '../modules/cutadapt/main.nf'
include { UMI_EXTRACT }               from '../modules/umi_extract/main.nf'
include { FASTQC }                    from '../modules/fastqc/main.nf'
include { BOWTIE2_FILTER }            from '../modules/bowtie2/main.nf'

workflow PREPROCESS_READS {

    take:
        // (sample_id, fastq, adapter_5, adapter_3, bc_pattern, has_umi)
        reads_ch
        contam_index

    main:
        /*
         * 0. Sanity check bc_pattern (only if has_umi)
         */
        checked = reads_ch.map { sid, fq, adapter_5, adapter_3, bc, has_umi ->
            if (has_umi) {
                assert bc && bc.contains('(?P<') && bc.startsWith('^') && bc.endsWith('$') :
                    "Malformed bc_pattern for sample ${sid}: '${bc}'"
            }
            tuple(sid, fq, adapter_5, adapter_3, bc, has_umi)
        }

        /*
         * 1. Adapter trimming
         */
        cutadapt_in = checked.map { sid, fq, adapter_5, adapter_3, bc, has_umi ->
            tuple(sid, fq, adapter_5, adapter_3, bc)
        }

        // CUTADAPT output: (sid, trimmed_fastq, bc_pattern)
        trimmed0 = CUTADAPT(cutadapt_in).trimmed_reads

        /*
         * 2. Re-attach has_umi flag by join on sample_id
         *    join → (sid, trimmed_fastq, bc_pattern, has_umi)
         */
        meta_has_umi = checked.map { sid, fq, adapter_5, adapter_3, bc, has_umi ->
            tuple(sid, has_umi)
        }
        trimmed = trimmed0.join(meta_has_umi, by: 0)

        /*
         * 3. Split by UMI presence
         */
        with_umi    = trimmed.filter  { sid, fq, bc, has_umi -> has_umi }
        without_umi = trimmed.filter  { sid, fq, bc, has_umi -> !has_umi }

        /*
         * 4a. UMI extraction (only if has_umi)
         */
        umi_res = UMI_EXTRACT(
            with_umi.map { sid, fq, bc, has_umi ->
                tuple(sid, fq, bc)
            }
        )
        umi_fastq_ch = umi_res.umi_fastq   // (sid, umi.fastq.gz)

        /*
         * 4b. Pass-through for non-UMI samples
         */
        no_umi_fastq_ch = without_umi.map { sid, fq, bc, has_umi ->
            tuple(sid, fq)
        }

        /*
         * 5. Pre-bowtie "processed" (still not contam-filtered yet)
         *    (sid, clean_fastq, has_umi)
         */
        processed0 = umi_fastq_ch
            .map { sid, fq -> tuple(sid, fq, true) }
            .mix(
                no_umi_fastq_ch.map { sid, fq -> tuple(sid, fq, false) }
            )

        /*
         * 6. Contaminant filtering (does NOT depend on junctions)
         *    BOWTIE2_FILTER input: (sid, fastq), idx_dir, idx_prefix
         *    output clean_reads: (sid, clean_fastq)
         */
        bowtie_clean = BOWTIE2_FILTER(
            processed0.map { sid, fq, has_umi -> tuple(sid, fq) },
            file(contam_index).parent,
            file(contam_index).getName()
        ).clean_reads

        /*
         * 7. Re-attach has_umi after bowtie
         */
        processed = bowtie_clean
            .join(processed0.map { sid, fq, has_umi -> tuple(sid, has_umi) }, by: 0)
            .map { sid, fq, has_umi -> tuple(sid, fq, has_umi) }

        /*
         * 8. QC (side-effect only)
         */
        FASTQC(processed.map { sid, fq, has_umi -> tuple(sid, fq) })

    emit:
        processed
}