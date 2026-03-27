nextflow.enable.dsl = 2

/*
 * ------------------------------------------------------------
 * Subworkflows / modules
 * ------------------------------------------------------------
 */
include { PREPROCESS_READS }   from '../subworkflows/preprocess_reads.nf'
include { ALIGN_RIBO_DEDUP }   from '../subworkflows/align_ribo_dedup.nf'
include { PARSE_RIBO_ALIGN_METRICS } from '../modules/parse_ribo_align_metrics/main.nf'
include { RIBO_QC_BASIC }      from '../subworkflows/ribo_qc_basic.nf'
include { CDS_QUANT }          from '../modules/cds_quant/main.nf'
include { PSITE_TRACK }        from '../modules/psite_track/main.nf'
include { TRANSLATED_ORFS }    from '../modules/translated_orfs/main.nf'

include { ALIGN_RNA_FASTP }    from '../subworkflows/align_rna_fastp.nf'
include { MERGE_JUNCTIONS }    from '../modules/merge_junctions/main.nf'
include { RNA_QUANT }          from '../modules/rna_quant/main.nf'
include { RNA_TRACK }          from '../modules/rna_track/main.nf'


/*
 * ------------------------------------------------------------
 * Sentinel file for "no sjdb"
 * ------------------------------------------------------------
 */
def NO_FILE = file("${projectDir}/assets/NO_FILE")


/*
 * ------------------------------------------------------------
 * Helper function to get organism-specific reference paths
 * ------------------------------------------------------------
 */
def getGenomeRefs(organism) {
    def refs = params.genomes?."${organism}"
    if (!refs) {
        log.warn "No genome configuration found for organism '${organism}'. Using default parameters."
        return [
            contam_index   : params.contam_index,
            star_index     : params.star_index,
            star_rna_index : params.star_rna_index,
            gtf            : params.gtf,
            genome_sizes   : params.genome_sizes
        ]
    }
    return refs
}


/*
 * ------------------------------------------------------------
 * Helper to build a safe sj key
 * ------------------------------------------------------------
 */
def makeSjKey(meta) {
    def org = (meta.organism ?: 'Unknown').toString().trim().replaceAll(/\s+/, '_')
    def sj  = (meta.sj_set_id ?: 'global').toString().trim().replaceAll(/\s+/, '_')
    return "${org}-${sj}".toString()
}


/*
 * ------------------------------------------------------------
 * Core ribo-only workflow
 * ------------------------------------------------------------
 */
workflow RiboSuite {

    take:
        reads_ch
        gtf

    main:
        processed = PREPROCESS_READS(reads_ch, params.contam_index)

        processed_with_sjdb = processed.map { sid, fq, has_umi ->
            tuple(sid, fq, has_umi, NO_FILE)
        }

        aligned = ALIGN_RIBO_DEDUP(
            processed_with_sjdb,
            params.star_index
        )

        qc = RIBO_QC_BASIC(aligned.genome_bam, gtf)

        psite_for_quant = qc.psite_offset_qc.map { sid, bam, bai, offsets ->
            tuple(sid, bam, bai, offsets)
        }

        cds_quant = CDS_QUANT(
            psite_for_quant.map { meta, bam, bai, offsets ->
                tuple(meta, bam, bai, meta.gtf, offsets)
            }
        )

        if (params.enable_translated_orfs) {
            translated_orfs = TRANSLATED_ORFS(
                qc.psite_offset_qc.map { sid, bam, bai, offsets -> tuple(sid, bam, bai, offsets) },
                gtf
            )
        }

        def do_psite_track = params.enable_psite_track && ((params.psite_track_mode ?: 'both') != 'none')
        if (do_psite_track) {
            psite_track = PSITE_TRACK(
                qc.psite_offset_qc.map { meta, bam, bai, offsets ->
                    tuple(meta, bam, offsets, meta.genome_sizes)
                }
            )
        }

    emit:
        genome_bam            = aligned.genome_bam
        rpf_length_qc         = qc.rpf_length_qc
        psite_offset_qc       = qc.psite_offset_qc
        frame_periodicity_tsv = qc.frame_periodicity_tsv
        frame_periodicity_png = qc.frame_periodicity_png
        metagene_qc           = qc.metagene_qc
        cds_quant

        translated_orfs_tsv = params.enable_translated_orfs ? translated_orfs.tsv : Channel.empty()
        psite_tracks        = (params.enable_psite_track && ((params.psite_track_mode ?: 'both') != 'none')) ? psite_track.psite_tracks_by_len : Channel.empty()
}


/*
 * ------------------------------------------------------------
 * Integrated RNA + Ribo workflow
 *   - RNA: fastp -> STAR -> SJ.out.tab
 *   - Merge SJ per (organism, sj_set_id) to STAR sjdb format
 *   - Ribo: attach sjdb by (organism, sj_set_id)
 * ------------------------------------------------------------
 */
workflow RiboSuiteIntegrated {

    take:
        samples_tsv
        gtf_unused

    main:

        /*
         * Parse unified samplesheet
         *
         * Expected columns:
         *   sample_id, fastq1, fastq2, adapter_5, adapter_3, umi_5, umi_3, assay, organism, sj_set_id
         */
        samples_ch = Channel
            .fromPath(samples_tsv)
            .splitCsv(header: true, sep: '\t')
            .map { row ->

                def adapter_5 = (row.adapter_5 && row.adapter_5 != '0' && row.adapter_5 != 'NA') ? row.adapter_5 : null
                def adapter_3 = (row.adapter_3 && row.adapter_3 != '0' && row.adapter_3 != 'NA') ? row.adapter_3 : null

                int umi5 = row.umi_5?.trim() ? row.umi_5.toInteger() : 0
                int umi3 = row.umi_3?.trim() ? row.umi_3.toInteger() : 0

                def bc_pattern = null
                if (umi5 > 0 || umi3 > 0) {
                    def parts = []
                    if (umi5 > 0) parts << "(?P<umi_1>.{${umi5}})"
                    parts << '.+'
                    if (umi3 > 0) parts << "(?P<umi_2>.{${umi3}})"
                    bc_pattern = '^' + parts.join('') + '$'
                }
                def has_umi = (bc_pattern != null)

                def assay = row.assay?.trim()?.toLowerCase()
                if (!assay) assay = 'ribo'

                def organism = row.organism?.toString()?.trim()
                if (!organism) organism = 'Unknown'

                def sj_set_id = row.sj_set_id?.toString()?.trim()
                if (!sj_set_id) sj_set_id = 'global'

                def refs = getGenomeRefs(organism)

                def fq1 = file(row.fastq1)
                def fq2 = (row.fastq2 && row.fastq2 != 'NA' && row.fastq2 != '')
                    ? file(row.fastq2)
                    : null

                if (assay == 'ribo' && fq2 != null) {
                    log.warn "Ribo sample ${row.sample_id}: fastq2 provided; ribo pipeline currently expects single-end reads (fastq1 will be used)."
                }

                def meta = [
                    sample_id     : row.sample_id.toString(),
                    organism      : organism,
                    sj_set_id     : sj_set_id,
                    is_paired     : (fq2 != null),
                    contam_index  : file(refs.contam_index),
                    star_index    : file(refs.star_index),
                    star_rna_index: file(refs.star_rna_index),
                    gtf           : file(refs.gtf),
                    genome_sizes  : file(refs.genome_sizes)
                ]

                tuple(
                    meta,
                    fq1,
                    fq2,
                    adapter_5,
                    adapter_3,
                    bc_pattern,
                    has_umi,
                    assay
                )
            }

        ribo_raw = samples_ch.filter { it[7] == 'ribo' }
        rna_raw  = samples_ch.filter { it[7] == 'rna'  }

        /*
         * Log organisms observed, but do not restrict to one per run
         */
        observed_orgs = samples_ch
            .map { meta, fq1, fq2, a5, a3, bc, has_umi, assay -> meta.organism }
            .distinct()
            .toList()

        observed_orgs.view { orgs ->
            "Detected organisms in run: ${orgs.join(', ')}"
        }

        /*
         * RNA: build meta + reads
         *   reads = fq1         for SE
         *   reads = [fq1, fq2]  for PE
         */
        rna_reads_ch = rna_raw.map { meta, fq1, fq2, a5, a3, bc, has_umi, assay ->
            def reads = meta.is_paired ? [fq1, fq2] : fq1
            tuple(meta, reads, meta.star_rna_index)
        }

        rna_aln = ALIGN_RNA_FASTP(rna_reads_ch)

        rna_bam_for_quant = rna_aln.bam.map { meta, bam ->
            tuple(meta, bam, meta.gtf)
        }

        rna_quant = RNA_QUANT(rna_bam_for_quant)

        def do_rna_track = ((params.rna_track_mode ?: 'both') != 'none')
        if (do_rna_track) {
            rna_track = RNA_TRACK(
                rna_aln.bam.map { meta, bam -> tuple(meta, bam, meta.genome_sizes) }
            )
        }

        /*
         * Merge RNA junctions by (organism, sj_set_id)
         */
        merged_sj = MERGE_JUNCTIONS(
            rna_aln.sj
                .map { meta, sj_tab ->
                    def sj_key = makeSjKey(meta)
                    tuple(sj_key, sj_tab)
                }
                .groupTuple(by: 0)
        )

        /*
         * Build default sjdb map for all ribo sj keys
         */
        sjdb_default_ch = ribo_raw
            .map { meta, fq1, fq2, a5, a3, bc, has_umi, assay ->
                makeSjKey(meta)
            }
            .distinct()
            .map { sj_key -> tuple(sj_key, NO_FILE) }

        sjdb_real_ch = merged_sj.sjdb
            .map { sj_key, sjdb_path -> tuple(sj_key.toString(), sjdb_path) }

        sjdb_by_set_ch = sjdb_default_ch
            .mix(sjdb_real_ch)
            .groupTuple(by: 0)
            .map { sj_key, paths ->
                def real = paths.find { it.getName() != 'NO_FILE' }
                tuple(sj_key, real ?: NO_FILE)
            }

        /*
         * Ribo preprocess
         */
        ribo_reads_ch = ribo_raw.map { meta, fq1, fq2, a5, a3, bc, has_umi, assay ->
            tuple(meta, fq1, a5, a3, bc, has_umi, meta.contam_index)
        }

        ribo_processed = PREPROCESS_READS(ribo_reads_ch)

        /*
         * Attach sjdb using organism-aware sj key
         */
        ribo_processed2 = ribo_processed.map { meta, clean_fq, has_umi ->
            def sj_key = makeSjKey(meta)
            tuple(meta, clean_fq, has_umi, sj_key)
        }

        sjdb_lookup = sjdb_by_set_ch
            .toList()
            .map { list ->
                def result = [:]
                list.each { item ->
                    result[item[0].toString()] = item[1]
                }
                result
            }

        ribo_for_align = ribo_processed2
            .combine(sjdb_lookup)
            .map { meta, clean_fq, has_umi, sj_key, sjdb_map ->
                def sjdb_path = sjdb_map[sj_key] ?: NO_FILE

                if (sjdb_path.getName() == 'NO_FILE') {
                    log.warn "Sample ${meta.sample_id}: no RNA junctions for sj_key='${sj_key}', aligning without sjdb"
                }

                tuple(meta, clean_fq, has_umi, sjdb_path, meta.star_index)
            }

        /*
         * Ribo alignment + dedup
         */
        aligned = ALIGN_RIBO_DEDUP(ribo_for_align)

        /*
         * Parse final STAR metrics
         */
        ribo_metrics_files = PARSE_RIBO_ALIGN_METRICS(aligned.align_log)

        ribo_metrics = ribo_metrics_files.metrics
            .splitCsv(header: true, sep: '\t')
            .filter { row ->
                row.unique_mapped_reads && row.unique_mapped_reads != 'NA'
            }
            .map { row ->
                tuple(
                    row.sample_id.toString(),
                    row.unique_mapped_reads.toString().trim().toLong(),
                    row.unique_mapped_pct.toString().trim().toBigDecimal()
                )
            }

        aligned_with_metrics = aligned.genome_bam
            .map { meta, bam, bai ->
                tuple(meta.sample_id.toString(), meta, bam, bai)
            }
            .join(ribo_metrics, by: 0)
            .map { sample_id, meta, bam, bai, unique_reads, unique_pct ->
                tuple(meta, bam, bai, unique_reads, unique_pct)
            }

        /*
         * Filter low-depth ribo samples
         */
        ribo_split = aligned_with_metrics.branch { meta, bam, bai, unique_reads, unique_pct ->
            pass: (!params.skip_low_depth_ribo) || (unique_reads >= params.min_unique_mapped_reads)
            fail: true
        }

        ribo_split.fail.view { meta, bam, bai, unique_reads, unique_pct ->
            "[SKIP] ${meta.sample_id}: uniquely mapped reads ${unique_reads} < ${params.min_unique_mapped_reads}, downstream ribo analysis skipped"
        }

        /*
         * Ribo QC
         */
        qc = RIBO_QC_BASIC(
            ribo_split.pass.map { meta, bam, bai, unique_reads, unique_pct ->
                tuple(meta, bam, bai, meta.gtf)
            }
        )

        /*
         * Filter P-site offset tables by number of lines
         * skip if lines <= params.min_psite_lines
         */
        /*
         * RIBO_QC_BASIC already filters invalid P-site offset tables internally
         */
        psite_for_quant = qc.psite_offset_qc.map { meta, bam, bai, offsets ->
            tuple(meta, bam, bai, offsets)
        }

        /*
         * Downstream ribo analysis only for samples passing both filters
         */
        cds_quant = CDS_QUANT(
            psite_for_quant.map { meta, bam, bai, offsets ->
                tuple(meta, bam, bai, meta.gtf, offsets)
            }
        )

        if (params.enable_translated_orfs) {
            translated_orfs = TRANSLATED_ORFS(
                psite_for_quant.map { meta, bam, bai, offsets ->
                    tuple(meta, bam, bai, offsets, meta.gtf)
                }
            )
        }

        def do_psite_track = params.enable_psite_track && ((params.psite_track_mode ?: 'both') != 'none')
        if (do_psite_track) {
            psite_track = PSITE_TRACK(
                psite_for_quant.map { meta, bam, bai, offsets ->
                    tuple(meta, bam, offsets, meta.genome_sizes)
                }
            )
        }

    emit:
        // RNA outputs
        rna_bam         = rna_aln.bam
        rna_sj          = rna_aln.sj
        rna_gene_counts = rna_aln.gene_counts
        rna_log         = rna_aln.log
        rna_fastp_json  = rna_aln.fastp_json
        rna_fastp_html  = rna_aln.fastp_html
        sjdb_by_set     = merged_sj.sjdb

        rna_fc_counts   = rna_quant.counts
        rna_fc_summary  = rna_quant.summary

        rna_track_all = do_rna_track ? rna_track.rna_track_all : Channel.empty()
        rna_track_pos = do_rna_track ? rna_track.rna_track_pos : Channel.empty()
        rna_track_neg = do_rna_track ? rna_track.rna_track_neg : Channel.empty()

        // Ribo outputs
        genome_bam            = aligned.genome_bam
        rpf_length_qc         = qc.rpf_length_qc
        psite_offset_qc       = psite_for_quant
        frame_periodicity_tsv = qc.frame_periodicity_tsv
        frame_periodicity_png = qc.frame_periodicity_png
        metagene_qc           = qc.metagene_qc
        cds_quant

        translated_orfs_tsv = params.enable_translated_orfs ? translated_orfs.tsv : Channel.empty()
        psite_tracks        = do_psite_track ? psite_track.psite_tracks_by_len : Channel.empty()

        low_depth_ribo = ribo_split.fail.map { meta, bam, bai, unique_reads, unique_pct ->
            tuple(meta, unique_reads, unique_pct)
        }
}


/*
 * ------------------------------------------------------------
 * Default entry workflow
 * ------------------------------------------------------------
 */
workflow {
    RiboSuiteIntegrated(
        params.samples,
        null
    )
}