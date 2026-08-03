process PROKKA_VERSION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Deliberately reusing the exact same container as modules/nf-core/prokka: its
    // environment.yml is kept in sync with the vendored module's (same prokka/openjdk/
    // parallel versions) specifically so this local version-reporting workaround can
    // share the vendored module's container instead of building/maintaining a near-
    // duplicate one.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7b/7bc89d4083c0a4baaaca0ef9ac0ba65e0feaebd8d88fe1c16c47041cbc67f360/data'
:         'community.wave.seqera.io/library/prokka_openjdk_parallel:f21b98bcef4c3579' }"

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
