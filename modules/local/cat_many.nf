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
    tuple val(meta), path("${prefix}.gz"), emit: concatenated_files
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
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
    pigz -c temp_concatenated > ${prefix}.gz

    echo "Cleaning up..."
    rm temp_concatenated

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gz
    echo "${task.process}:" > versions.yml
    echo "    pigz: 2.6" >> versions.yml
    echo "    gzip: 1.10" >> versions.yml
    """
}
