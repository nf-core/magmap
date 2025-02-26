process RENAME_CONTIGS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_renamed_contigs.fna.gz"), emit: renamed_contigs
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id

    """
    seqkit replace -p "^" -r "${prefix}_" $contigs | gzip -c > ${prefix}_renamed_contigs.fna.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//' | sed 's/ Build.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_renamed_contigs.fna.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.3.1
    END_VERSIONS
    """
}
