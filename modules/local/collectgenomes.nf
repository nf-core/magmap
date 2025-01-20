process COLLECTGENOMES {
    tag "$accno"
    label 'process_single'

    conda "conda-forge::wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hed695b0_4':
        'biocontainers/gnu-wget:1.18--hed695b0_4' }"

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
        link_2="\$(grep "$accno" $genomeinfo | cut -d "," -f 2)"
        link_3="\$(grep "$accno" $genomeinfo | cut -d "," -f 3)"

        # Check whether link_2 is a web link ("http" or "https")
        if [[ \$link_2 == http* || \$link_2 == https* ]]; then
            echo "Downloading \$link_2"
            wget "\$link_2"
        else
            echo "Creating symlink for \$link_2"
            ln -s "\$link_2" .
        fi

        # Check whether link_3 is a web link ("http" or "https")
        if [[ \$link_3 == http* || \$link_3 == https* ]]; then
            echo "Downloading \$link_3"
            wget "\$link_3"
        else
            echo "Creating symlink for \$link_3"
            ln -s "\$link_3" .
        fi
    else
        echo "The string is not present in the file."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | grep 'GNU Wget' | sed 's/GNU Wget \\([0-9.]\\+\\) .*/\\1/')
    END_VERSIONS
    """
}
