process KOFAMSCAN_TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(tsvs, stageAs: 'tsvs/')

    output:
    tuple val(meta), path("${prefix}.kofamscan.tsv.gz"), emit: tsv
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "orf\tko\tthrshld\tscore\tevalue\tko_definition" | gzip -c > ${prefix}.kofamscan.tsv.gz

    find tsvs -name "*.tsv.gz" -print0 \\
        | sort -z \\
        | xargs -0 zcat \\
        | { grep -v '^#' || true; } \\
        | cut -f 2-7 \\
        | sed 's/\\t"/\\t/' \\
        | sed 's/"\$//' \\
        | gzip -c >> ${prefix}.kofamscan.tsv.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "orf\tko\tthrshld\tscore\tevalue\tko_definition" | gzip -c > ${prefix}.kofamscan.tsv.gz
    """
}
