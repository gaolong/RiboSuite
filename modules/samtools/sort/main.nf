nextflow.enable.dsl=2

process SAMTOOLS_SORT {

    tag "$sample_id"

    cpus   { params.samtools_sort_cpus ?: 1 }
    memory { params.samtools_sort_memory ?: '32 GB' }

    conda "bioconda::samtools=1.18"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id),
              path("${sample_id}.sorted.bam")

    script:
        """
        samtools sort \
          -@ ${task.cpus} \
          -o ${sample_id}.sorted.bam \
          $bam
        """
}
