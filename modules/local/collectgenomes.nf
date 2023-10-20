process COLLECTGENOMES {
    tag "$accno"
    label 'process_single'

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0':
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' }"

    input:
    val accno
    path genomeinfo

    output:
    tuple val(accno), path("*.fna.gz"), optional: true, emit: fna
    tuple val(accno), path("*.gff.gz"), optional: true, emit: gff
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    if grep -q $accno $genomeinfo; then
        echo "The string is present in the file."
        ln -s "\$(grep "$accno" $genomeinfo | cut -d "," -f 2)" .
        ln -s "\$(grep "$accno" $genomeinfo | cut -d "," -f 3)" .
    else
        echo "The string is not present in the file."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
