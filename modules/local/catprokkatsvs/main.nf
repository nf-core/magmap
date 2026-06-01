process CATPROKKATSVS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8':
        'biocontainers/pigz:2.8' }"

    input:
    tuple val(meta), path(tsvs, stageAs: 'tsvs/')

    output:
    tuple val(meta), path("${outfile}.gz"), emit: tsv
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^pigz //"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outfile = "${prefix}.prokka-annotations.tsv"
    """
    echo "orf\tftype\tlength\tgene\tec_number\tcog\tproduct" > ${outfile}

    find tsvs -name "*.tsv.gz" | xargs unpigz -c | grep -v '^locus_tag' | grep -v '\tgene\t' >> ${outfile}

    pigz ${outfile}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    outfile = "${prefix}.prokka-annotations.tsv"
    """
    echo $args
    touch ${outfile}.gz
    """
}
