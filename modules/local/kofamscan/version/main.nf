process KOFAMSCAN_VERSION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Keep environment.yml identical to modules/nf-core/kofamscan so both share one container.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2':
        'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

    output:
    tuple val("${task.process}"), val('kofamscan'), eval("exec_annotation --version 2>&1 | sed 's/exec_annotation //;'"), emit: versions_kofamscan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """

    stub:
    """
    """
}
