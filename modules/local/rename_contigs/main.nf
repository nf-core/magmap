process RENAME_CONTIGS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${prefix}.renamed.fna.gz"), emit: renamed_contigs
    tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//" | sed "s/ Build.*//"'), emit: versions_seqkit, topic: versions

    script:
    prefix     = task.ext.prefix ?: meta.id
    def prefix_md5 = prefix.md5().substring(0,9)

    """
    seqkit replace -p "^" -r "${prefix_md5}_" $contigs | gzip -c > ${prefix}.renamed.fna.gz
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    cat /dev/null | gzip -c > ${prefix}.renamed.fna.gz
    """
}
