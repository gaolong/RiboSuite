nextflow.enable.dsl=2

process STAR_RIBO_ALIGN {

    tag "$sample_id"

    conda "bioconda::star=2.7.11b"

    input:
        tuple val(sample_id), path(reads)
        path star_index

    output:
        tuple val(sample_id),
              path("${sample_id}.Aligned.sortedByCoord.out.bam"),
              path("${sample_id}.Aligned.toTranscriptome.out.bam"), optional: true

    script:
        def quantArg = params.star_quantmode ? "--quantMode TranscriptomeSAM" : ""

        """
        STAR \
          --genomeDir ${star_index} \
          --readFilesIn ${reads} \
          --readFilesCommand zcat \
          --runThreadN ${task.cpus} \
          --outSAMtype BAM SortedByCoordinate \
          ${quantArg} \
          --outFileNamePrefix ${sample_id}.
        """
}
