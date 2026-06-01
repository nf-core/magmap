process CAT_MANY {
    tag "${meta.id}"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' }"

    input:
    val(meta)
    path(files2cat), stageAs: 'input/*'

    output:
    tuple val(meta), path("*.gz"), emit: concatenated_files
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^pigz //"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def outfile = "${prefix}.gz"

    """
    echo "Concatenating files..."
    for file in ./input/*; do
        if [ -f "\$file" ]; then
            if [[ "\$file" == *.gz ]]; then
                zcat "\$file" >> temp_concatenated
            else
                cat "\$file" >> temp_concatenated
            fi
        fi
    done

    echo "Compressing concatenated file..."
    pigz -c temp_concatenated > ${outfile}

    echo "Cleaning up..."
    rm temp_concatenated
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gz
    """
}
