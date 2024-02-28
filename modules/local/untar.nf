process UNTAR {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::aria2=1.36.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aria2:1.36.0' :
        'biocontainers/aria2:1.36.0' } "

    input:
    tuple val(meta), path(tar)

    output:
    tuple val(meta), path ("checkm_data_2015_01_16/"), emit: downloaded_file
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir checkm_data_2015_01_16/
    tar x -C checkm_data_2015_01_16 -v -z -f *.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | sed 's/bsdtar //g' | sed 's/- libarchive 3.5.3 zlib\\/1.2.12 liblzma\\/5.0.5 bz2lib\\/1.0.8 // )
    END_VERSIONS
    """
}
