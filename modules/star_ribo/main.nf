process STAR_RIBO_ALIGN {

    tag "$sample_id"

    conda "bioconda::star=2.7.11b"

    publishDir "${params.outdir}/align/star_ribo",
               mode: 'copy',
               saveAs: { file -> "${sample_id}/${file}" }

    input:
        tuple val(sample_id), path(reads)
        path star_index

    output:
        /*
         * Canonical genome alignment (always present)
         */
        tuple val(sample_id),
              path("${sample_id}.Aligned.sortedByCoord.out.bam"),
              emit: genome_bam

        /*
         * Transcriptome alignment (only when --star_quantmode true)
         */
        tuple val(sample_id),
              path("${sample_id}.Aligned.toTranscriptome.out.bam"),
              optional: true,
              emit: tx_bam

        /*
         * STAR splice junctions (always produced)
         */
        path "${sample_id}.SJ.out.tab", emit: sj_tab

        /*
         * STAR final log (always useful)
         */
        path "${sample_id}.Log.final.out", emit: star_log

    script:
        def quantArg = params.star_quantmode ? "--quantMode TranscriptomeSAM" : ""

        """
        set -euo pipefail

        STAR \
            --genomeDir ${star_index} \
            --readFilesIn ${reads} \
            --readFilesCommand zcat \
            --runThreadN ${task.cpus} \
            --alignEndsType EndToEnd \
            --outFilterMismatchNmax 2 \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 51 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All \
            ${quantArg} \
            --outTmpDir STARtmp \
            --outFileNamePrefix ${sample_id}.
        """
}
