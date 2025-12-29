nextflow.enable.dsl=2

process STAR_ALIGN {

    tag "$sample_id"

    cpus   { params.star_cpus ?: 1 }
    memory { params.star_memory ?: '32 GB' }

    conda "bioconda::star=2.7.11b"

    input:
        tuple val(sample_id), path(reads)
        path star_index

    output:
        tuple val(sample_id),
              path("${sample_id}.Aligned.sortedByCoord.out.bam")

    script:
        """
        STAR \
          --genomeDir ${star_index} \
          --readFilesIn ${reads} \
          --readFilesCommand zcat \
          --runThreadN ${task.cpus} \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix ${sample_id}.
        """
}
