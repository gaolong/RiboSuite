nextflow.enable.dsl = 2

process STAR_RIBO_ALIGN {

    tag "${meta.sample_id}"

    conda "bioconda::star=2.7.11b"

    publishDir "${params.outdir}/ribo/align/star_ribo",
        mode: 'copy',
        saveAs: { file -> "${meta.sample_id}/${file}" }

    input:
        tuple val(meta), path(reads), val(has_umi), path(sjdb), path(star_index)

    output:
        tuple val(meta),
            path("${meta.sample_id}.Aligned.sortedByCoord.out.bam"),
            emit: genome_bam

        tuple val(meta),
            path("${meta.sample_id}.Aligned.toTranscriptome.out.bam"),
            optional: true,
            emit: tx_bam

        tuple val(meta),
            path("${meta.sample_id}.SJ.out.tab"),
            emit: sj_tab

        tuple val(meta),
            path("${meta.sample_id}.Log.final.out"),
            emit: star_log

    script:

        def quantModeOn = (params.star_quantmode ?: false) as boolean
        def quantArg    = quantModeOn ? "--quantMode TranscriptomeSAM" : ""

        def rescue_thr  = (params.star_rescue_unique_pct ?: 30)

        def sjdbProvided = (sjdb && sjdb.getName() != 'NO_FILE')
        def sjdbArg      = sjdbProvided ? "--sjdbFileChrStartEnd ${sjdb}" : ""

        def readsArg = (reads instanceof List)
            ? reads.join(' ')
            : reads.toString()

        """
        set -euo pipefail

        ########################################
        # 1. Primary alignment: EndToEnd
        ########################################
        STAR \\
            --genomeDir ${star_index} \\
            --readFilesIn ${readsArg} \\
            --readFilesCommand zcat \\
            --runThreadN ${task.cpus} \\
            ${sjdbArg} \\
            --alignEndsType EndToEnd \\
            --outFilterMismatchNmax 2 \\
            --alignSJDBoverhangMin 1 \\
            --alignSJoverhangMin 51 \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMattributes All \\
            ${quantArg} \\
            --outTmpDir STARtmp.EndToEnd \\
            --outFileNamePrefix ${meta.sample_id}.

        ########################################
        # 2. Parse uniquely mapped rate
        ########################################
        uniq_pct=\$(grep -F "Uniquely mapped reads %" ${meta.sample_id}.Log.final.out | awk '{print \$6}')
        uniq_pct=\${uniq_pct%\\%}

        echo "[STAR_RIBO_ALIGN] ${meta.sample_id}: uniquely mapped = \${uniq_pct}%" >> ${meta.sample_id}.Log.final.out
        echo "[STAR_RIBO_ALIGN] ${meta.sample_id}: sjdb injected = ${sjdbProvided ? "YES" : "NO"}" >> ${meta.sample_id}.Log.final.out
        echo "[STAR_RIBO_ALIGN] ${meta.sample_id}: quantMode TranscriptomeSAM = ${quantModeOn ? "ON" : "OFF"}" >> ${meta.sample_id}.Log.final.out

        ########################################
        # 3. Rescue if needed
        ########################################
        if awk -v a="\$uniq_pct" -v b="${rescue_thr}" 'BEGIN{exit !(a < b)}'; then

            echo "" >> ${meta.sample_id}.Log.final.out
            echo "===== RESCUE ALIGNMENT TRIGGERED =====" >> ${meta.sample_id}.Log.final.out
            echo "Original alignEndsType: EndToEnd" >> ${meta.sample_id}.Log.final.out
            echo "Rescue alignEndsType:   Local" >> ${meta.sample_id}.Log.final.out
            echo "Uniquely mapped reads:  \${uniq_pct}%" >> ${meta.sample_id}.Log.final.out
            echo "STAR sjdb injected:     ${sjdbProvided ? "YES" : "NO"}" >> ${meta.sample_id}.Log.final.out
            echo "STAR quantMode:         ${quantModeOn ? "TranscriptomeSAM" : "OFF"}" >> ${meta.sample_id}.Log.final.out
            echo "=====================================" >> ${meta.sample_id}.Log.final.out
            echo "" >> ${meta.sample_id}.Log.final.out

            STAR \\
                --genomeDir ${star_index} \\
                --readFilesIn ${readsArg} \\
                --readFilesCommand zcat \\
                --runThreadN ${task.cpus} \\
                ${sjdbArg} \\
                --alignEndsType Local \\
                --outFilterMismatchNmax 2 \\
                --alignSJDBoverhangMin 1 \\
                --alignSJoverhangMin 51 \\
                --outSAMtype BAM SortedByCoordinate \\
                --outSAMattributes All \\
                ${quantArg} \\
                --outTmpDir STARtmp.Local \\
                --outFileNamePrefix ${meta.sample_id}.rescue.

            mv ${meta.sample_id}.rescue.Aligned.sortedByCoord.out.bam \\
               ${meta.sample_id}.Aligned.sortedByCoord.out.bam

            mv ${meta.sample_id}.rescue.SJ.out.tab \\
               ${meta.sample_id}.SJ.out.tab

            if [ -f ${meta.sample_id}.rescue.Aligned.toTranscriptome.out.bam ]; then
                mv ${meta.sample_id}.rescue.Aligned.toTranscriptome.out.bam \\
                   ${meta.sample_id}.Aligned.toTranscriptome.out.bam
            fi

            echo "" >> ${meta.sample_id}.Log.final.out
            echo "===== RESCUE STAR Log.final.out =====" >> ${meta.sample_id}.Log.final.out
            cat ${meta.sample_id}.rescue.Log.final.out >> ${meta.sample_id}.Log.final.out
        fi
        """
}