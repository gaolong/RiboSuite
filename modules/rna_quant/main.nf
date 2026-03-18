nextflow.enable.dsl = 2

process RNA_QUANT {

    tag "${meta?.sample_id ?: meta?.id ?: 'sample'}"

    conda "bioconda::subread=2.0.6"

    publishDir "${params.outdir}/rna/quant/",
        mode: 'copy',
        saveAs: { file -> "${meta?.sample_id ?: meta?.id ?: 'sample'}/${file}" }

    input:
        tuple val(meta), path(bam), path(gtf)

    output:
        tuple val(meta), path("${meta?.sample_id ?: meta?.id ?: 'sample'}.featureCounts.txt"), emit: counts
        tuple val(meta), path("${meta?.sample_id ?: meta?.id ?: 'sample'}.featureCounts.txt.summary"), emit: summary

    script:
        def sample_id    = meta?.sample_id ?: meta?.id ?: 'sample'
        def stranded     = (params.rna_quant_strandedness != null) ? params.rna_quant_strandedness as Integer : 0
        def is_paired    = meta?.containsKey('is_paired') ? meta.is_paired : false
        def paired       = is_paired ? "-p" : ""
        def feature_type = params.rna_quant_feature_type ?: "exon"
        def gene_attr    = params.rna_quant_gene_attr ?: "gene_id"
        def extra_args   = params.rna_quant_extra_args ?: ""

        """
        featureCounts \\
            -T ${task.cpus} \\
            -a ${gtf} \\
            -o ${sample_id}.featureCounts.txt \\
            -t ${feature_type} \\
            -g ${gene_attr} \\
            -s ${stranded} \\
            ${paired} \\
            ${extra_args} \\
            ${bam}
        """
}