process CAT_MANY {
    label 'process_long'

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' }"

    input:
    val(meta)
    path(files2cat), stageAs: 'input/*'

    output:
    tuple val(meta), path("${meta.id}.gz"), emit: concatenated_files
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args   ?: ''
    def checkgz = file2cat.toString() =~ /\.gz$/ ? "zcat" : "cat"

    """
    find ./input/* -name '*' | xargs ${checkgz} | gzip -c >> ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
