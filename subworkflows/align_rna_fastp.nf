nextflow.enable.dsl = 2

include { FASTP }          from '../modules/fastp/main.nf'
include { STAR_RNA_ALIGN } from '../modules/star_rna/main.nf'

workflow ALIGN_RNA_FASTP {

    take:
        // tuple(meta, reads, star_index)
        reads_ch

    main:
        /*
         * FASTP only needs meta + reads
         */
        fp_input = reads_ch.map { meta, reads, star_index ->
            tuple(meta, reads)
        }

        fp = FASTP(fp_input)

        /*
         * Keep a lookup from sample_id -> star_index
         * so we can reattach per-sample STAR index after FASTP
         */
        star_index_by_sample = reads_ch.map { meta, reads, star_index ->
            tuple(meta.sample_id.toString(), star_index)
        }

        /*
         * Normalize FASTP output so reads is always a list
         */
        fp_reads_norm = fp.reads.map { meta, reads ->
            def read_list = (reads instanceof List) ? reads : [reads]
            tuple(meta, read_list)
        }

        /*
         * Reattach star_index by sample_id
         * Output shape: tuple(meta, cleaned_reads, star_index)
         */
        fp_reads_for_star = fp_reads_norm
            .map { meta, reads ->
                tuple(meta.sample_id.toString(), meta, reads)
            }
            .join(star_index_by_sample, by: 0)
            .map { sample_id, meta, reads, star_index ->
                tuple(meta, reads, star_index)
            }

        aln = STAR_RNA_ALIGN(fp_reads_for_star)

        /*
         * Normalize outputs to plain Nextflow tuples
         */
        bam_ch = aln.bam.map { meta, bam ->
            tuple(meta, bam)
        }

        sj_ch = aln.sj.map { meta, sj ->
            tuple(meta, sj)
        }

        gc_ch = aln.gene_counts.map { meta, gc ->
            tuple(meta, gc)
        }

        log_ch = aln.log.map { meta, log ->
            tuple(meta, log)
        }

    emit:
        bam           = bam_ch
        sj            = sj_ch
        gene_counts   = gc_ch
        star_log      = log_ch
        fastp_json    = fp.json
        fastp_html    = fp.html
        cleaned_reads = fp_reads_norm
}