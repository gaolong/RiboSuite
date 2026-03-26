nextflow.enable.dsl = 2

include { PSITE_SHIFT_BAM } from '../modules/psite_shift_bam/main.nf'
include { PSITE_TRACK     } from '../modules/psite_track/main.nf'
include { CDS_QUANT       } from '../modules/cds_quant/main.nf'
include { TRANSLATED_ORFS } from '../modules/translated_orfs/main.nf'

params.samples          = params.samples ?: null
params.genome           = params.genome ?: null
params.outdir           = params.outdir ?: "results_pool_from_bam"
params.psite_track_mode = params.psite_track_mode ?: "both"
params.psite_by_length  = params.psite_by_length ?: false

process MERGE_BAM {

    tag "${meta.merge_id}"

    conda "bioconda::samtools=1.18"

    publishDir "${params.outdir}/ribo/pool/bam",
        mode: 'copy',
        saveAs: { file -> "${meta.merge_id}/${file}" }

    input:
        tuple val(meta), path(bams)

    output:
        tuple val(meta),
              path("${meta.merge_id}.merged.bam"),
              path("${meta.merge_id}.merged.bam.bai")

    script:
    def bam_args = bams.collect { it.getName() }.join(' ')

    """
    set -euo pipefail

    samtools merge -f ${meta.merge_id}.merged.bam ${bam_args}
    samtools index ${meta.merge_id}.merged.bam
    """
}

process MAKE_ZERO_OFFSET_TABLE {

    tag "${meta.merge_id}"

    conda "conda-forge::python=3.11"

    publishDir "${params.outdir}/ribo/pool/offset",
        mode: 'copy',
        saveAs: { file -> "${meta.merge_id}/${file}" }

    input:
        tuple val(meta), path(offset_tables)

    output:
        tuple val(meta), path("${meta.merge_id}.psite_shifted.offsets.tsv")

    script:
    """
    set -euo pipefail

    python - << 'PY'
import csv

files = "${offset_tables.join(' ')}".split()
lengths = set()

for f in files:
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        if 'read_length' not in reader.fieldnames:
            raise ValueError(f"Missing read_length column in {f}")
        for row in reader:
            lengths.add(int(row['read_length']))

lengths = sorted(lengths)

with open("${meta.merge_id}.psite_shifted.offsets.tsv", "w", newline="") as fh:
    writer = csv.writer(fh, delimiter='\\t')
    writer.writerow(["sample_id", "read_length", "psite_offset"])
    for rl in lengths:
        writer.writerow(["${meta.merge_id}", rl, 0])
PY
    """
}

workflow {

    if( !params.samples ) {
        error "Missing required parameter: --samples"
    }

    if( !params.genome ) {
        error "Missing required parameter: --genome"
    }

    if( !params.genomes ) {
        error "Missing params.genomes in nextflow.config"
    }

    if( !params.genomes.containsKey(params.genome) ) {
        error "Genome '${params.genome}' not found in params.genomes"
    }

    def genome_cfg   = params.genomes[params.genome]
    def genome_sizes = genome_cfg?.genome_sizes
    def gtf          = genome_cfg?.gtf

    if( !genome_sizes ) {
        error "Missing genome_sizes for genome '${params.genome}' in nextflow.config"
    }

    if( !gtf ) {
        error "Missing gtf for genome '${params.genome}' in nextflow.config"
    }

    /*
     * TSV columns:
     * sample_id  merge_id  bam  bai  offset
     */
    ch_rows_all = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            assert row.sample_id : "Missing sample_id: ${row}"
            assert row.merge_id  : "Missing merge_id: ${row}"
            assert row.bam       : "Missing bam: ${row}"
            assert row.bai       : "Missing bai: ${row}"
            assert row.offset    : "Missing offset: ${row}"

            def meta = [
                sample_id: row.sample_id.toString(),
                merge_id : row.merge_id.toString()
            ]

            tuple(
                meta,
                file(row.bam,    checkIfExists: true),
                file(row.bai,    checkIfExists: true),
                file(row.offset, checkIfExists: true)
            )
        }

    /*
     * ------------------------------------------------------------
     * Step 1: unique samples only for PSITE_SHIFT_BAM
     * ------------------------------------------------------------
     */
    ch_unique_samples = ch_rows_all
        .map { meta, bam, bai, offset ->
            tuple(meta.sample_id, bam, bai, offset)
        }
        .unique { it[0] }
        .map { sample_id, bam, bai, offset ->
            tuple(
                [sample_id: sample_id],
                bam,
                bai,
                offset
            )
        }

    PSITE_SHIFT_BAM(ch_unique_samples)

    /*
     * ------------------------------------------------------------
     * Step 2: Create a map for easy lookup
     * ------------------------------------------------------------
     * Instead of join, we'll use combine and filter
     */
    ch_membership = ch_rows_all
        .map { meta, bam, bai, offset ->
            tuple(meta.sample_id, meta.merge_id)
        }

    ch_shifted = PSITE_SHIFT_BAM.out.shifted
        .map { meta, shifted_bam, shifted_bai, shifted_offsets ->
            tuple(meta.sample_id, shifted_bam, shifted_bai, shifted_offsets)
        }

    /*
     * Use combine instead of join to preserve all membership relationships
     */
    ch_membership_with_shifted = ch_membership
        .combine(ch_shifted, by: 0)  // Combine by sample_id (index 0)
        .map { sample_id, merge_id, shifted_bam, shifted_bai, shifted_offsets ->
            tuple(merge_id, shifted_bam, shifted_offsets)
        }

    /*
     * ------------------------------------------------------------
     * Step 3: Split into separate channels using multiMap
     * ------------------------------------------------------------
     */
    ch_split = ch_membership_with_shifted
        .multiMap { merge_id, shifted_bam, shifted_offsets ->
            bams: tuple(merge_id, shifted_bam)
            offsets: tuple(merge_id, shifted_offsets)
        }

    /*
     * ------------------------------------------------------------
     * Step 4: Group by merge_id
     * ------------------------------------------------------------
     */
    ch_shifted_bams_grouped = ch_split.bams
        .groupTuple()
        .map { merge_id, bam_list ->
            tuple([merge_id: merge_id, sample_id: merge_id], bam_list)
        }

    ch_shifted_offsets_grouped = ch_split.offsets
        .groupTuple()
        .map { merge_id, offset_list ->
            tuple([merge_id: merge_id, sample_id: merge_id], offset_list)
        }

    MERGE_BAM(ch_shifted_bams_grouped)
    MAKE_ZERO_OFFSET_TABLE(ch_shifted_offsets_grouped)

    /*
     * ------------------------------------------------------------
     * Step 5: join merged shifted BAM with merged zero-offset table
     * ------------------------------------------------------------
     */
    ch_merged_joined = MERGE_BAM.out
        .map { meta, merged_bam, merged_bai ->
            tuple(meta.merge_id, meta, merged_bam, merged_bai)
        }
        .join(
            MAKE_ZERO_OFFSET_TABLE.out.map { meta, zero_offsets ->
                tuple(meta.merge_id, zero_offsets)
            }
        )
        .map { merge_id, meta, merged_bam, merged_bai, zero_offsets ->
            tuple(meta, merged_bam, merged_bai, zero_offsets)
        }

    /*
     * ------------------------------------------------------------
     * Step 6: pooled PSITE_TRACK
     * ------------------------------------------------------------
     */
    ch_psite_input = ch_merged_joined
        .map { meta, merged_bam, merged_bai, zero_offsets ->
            tuple(
                meta,
                merged_bam,
                zero_offsets,
                file(genome_sizes, checkIfExists: true)
            )
        }

    PSITE_TRACK(ch_psite_input)

    /*
     * ------------------------------------------------------------
     * Step 7: pooled CDS_QUANT
     * ------------------------------------------------------------
     */
    ch_cds_quant_input = ch_merged_joined
        .map { meta, merged_bam, merged_bai, zero_offsets ->
            tuple(
                meta,
                merged_bam,
                merged_bai,
                file(gtf, checkIfExists: true),
                zero_offsets
            )
        }

    CDS_QUANT(ch_cds_quant_input)

    /*
     * ------------------------------------------------------------
     * Step 8: pooled TRANSLATED_ORFS
     * ------------------------------------------------------------
     */
    ch_translated_orfs_input = ch_merged_joined
        .map { meta, merged_bam, merged_bai, zero_offsets ->
            tuple(
                meta,
                merged_bam,
                merged_bai,
                zero_offsets,
                file(gtf, checkIfExists: true)
            )
        }

    TRANSLATED_ORFS(ch_translated_orfs_input)
}