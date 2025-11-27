process CHECK_DUPLICATES {
    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'biocontainers/ubuntu:24.04' }"

    input:
    tuple val(meta), path(fnas, stageAs: 'contigs/*')

    output:
    stdout                                emit: duplicate_genomes
    path "*.genomes_with_duplicates.txt", emit: genomes_with_duplicates
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: meta.id
    def outfile = "${prefix}.genomes_with_duplicates.txt"

    """
    zgrep -h '>' $fnas | sed 's/>//' | sort | uniq -d > dupl_contig_names.txt
    zgrep -l -F -f dupl_contig_names.txt $fnas | sed 's:contigs/::' | sort -u > ${outfile} || touch ${outfile}
    cat ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zgrep: \$( zgrep --version | sed 's/.*/1.5/' | head -n 1 )
    END_VERSIONS
    """
}
