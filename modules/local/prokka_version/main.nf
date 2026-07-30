process PROKKA_VERSION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3af46b047c8fe84112adeaecf300878217c629b97f111f923ecf327656ddd141/data' :
        'community.wave.seqera.io/library/prokka_openjdk:10546cadeef11472' }"

    output:
    tuple val("${task.process}"), val('prokka'), eval("prokka --version 2>&1 | sed 's/^.*prokka //'"), emit: versions_prokka, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """

    stub:
    """
    """
}
