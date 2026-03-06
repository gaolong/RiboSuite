process STAR_RIBO_ALIGN {

    tag "${sample_id}"

    conda "bioconda::star=2.7.11b"

    publishDir "${params.outdir}/ribo/align/star_ribo",
        mode: 'copy',
        saveAs: { file -> "${sample_id}/${file}" }

    input:
        tuple val(sample_id), path(reads)
        path star_index
        path sjdb

    output:
        tuple val(sample_id),
            path("${sample_id}.Aligned.sortedByCoord.out.bam"),
            emit: genome_bam

        tuple val(sample_id),
            path("${sample_id}.Aligned.toTranscriptome.out.bam"),
            optional: true,
            emit: tx_bam

        path "${sample_id}.SJ.out.tab", emit: sj_tab
        path "${sample_id}.Log.final.out", emit: star_log

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
            --outFileNamePrefix ${sample_id}.

        ########################################
        # 2. Parse uniquely mapped rate
        ########################################
        uniq_pct=\$(grep -F "Uniquely mapped reads %" ${sample_id}.Log.final.out | awk '{print \$6}')
        uniq_pct=\${uniq_pct%\\%}

        echo "[STAR_RIBO_ALIGN] ${sample_id}: uniquely mapped = \${uniq_pct}%" >> ${sample_id}.Log.final.out
        echo "[STAR_RIBO_ALIGN] ${sample_id}: sjdb injected = ${sjdbProvided ? "YES" : "NO"}" >> ${sample_id}.Log.final.out
        echo "[STAR_RIBO_ALIGN] ${sample_id}: quantMode TranscriptomeSAM = ${quantModeOn ? "ON" : "OFF"}" >> ${sample_id}.Log.final.out

        ########################################
        # 3. Rescue if needed
        ########################################
        if awk -v a="\$uniq_pct" -v b="${rescue_thr}" 'BEGIN{exit !(a < b)}'; then

            echo "" >> ${sample_id}.Log.final.out
            echo "===== RESCUE ALIGNMENT TRIGGERED =====" >> ${sample_id}.Log.final.out
            echo "Original alignEndsType: EndToEnd" >> ${sample_id}.Log.final.out
            echo "Rescue alignEndsType:   Local" >> ${sample_id}.Log.final.out
            echo "Uniquely mapped reads:  \${uniq_pct}%" >> ${sample_id}.Log.final.out
            echo "STAR sjdb injected:     ${sjdbProvided ? "YES" : "NO"}" >> ${sample_id}.Log.final.out
            echo "STAR quantMode:         ${quantModeOn ? "TranscriptomeSAM" : "OFF"}" >> ${sample_id}.Log.final.out
            echo "=====================================" >> ${sample_id}.Log.final.out
            echo "" >> ${sample_id}.Log.final.out

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
                --outFileNamePrefix ${sample_id}.rescue.

            mv ${sample_id}.rescue.Aligned.sortedByCoord.out.bam \\
               ${sample_id}.Aligned.sortedByCoord.out.bam

            mv ${sample_id}.rescue.SJ.out.tab \\
               ${sample_id}.SJ.out.tab

            if [ -f ${sample_id}.rescue.Aligned.toTranscriptome.out.bam ]; then
                mv ${sample_id}.rescue.Aligned.toTranscriptome.out.bam \\
                   ${sample_id}.Aligned.toTranscriptome.out.bam
            fi

            echo "" >> ${sample_id}.Log.final.out
            echo "===== RESCUE STAR Log.final.out =====" >> ${sample_id}.Log.final.out
            cat ${sample_id}.rescue.Log.final.out >> ${sample_id}.Log.final.out
        fi
        """
}