process COLLECTGENOMES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0':
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' }"

    input:
    tuple val(meta), val(accno)
    tuple val(meta), path(fnas), path(gffs)

    output:
    tuple val(meta), path("*.fna.gz") , emit: fnas
    tuple val(meta), path("*.gff.gz") , emit: gffs
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Set the target filename to search for
    target_filename1="${accno}*.fna.gz"

    # Use the find command to search for the target filename
    found_file=\$(find . -name "$target_filename" -print -quit)

    # Check if the file was found
    mkdir genome
    if [[ -n "$found_file" ]]; then
        echo "File found: $found_file"
    # Create a symlink for the found file
        ln -s "$found_file" ./genome
    else
        echo "No matching file found"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
