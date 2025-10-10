process PROKKA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3af46b047c8fe84112adeaecf300878217c629b97f111f923ecf327656ddd141/data' :
        'community.wave.seqera.io/library/prokka_openjdk:10546cadeef11472' }"

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}/${prefix}.gff.gz"), emit: gff
    tuple val(meta), path("${prefix}/${prefix}.gbk.gz"), emit: gbk
    tuple val(meta), path("${prefix}/${prefix}.fna.gz"), emit: fna
    tuple val(meta), path("${prefix}/${prefix}.faa.gz"), emit: faa
    tuple val(meta), path("${prefix}/${prefix}.ffn.gz"), emit: ffn
    tuple val(meta), path("${prefix}/${prefix}.sqn.gz"), emit: sqn
    tuple val(meta), path("${prefix}/${prefix}.fsa.gz"), emit: fsa
    tuple val(meta), path("${prefix}/${prefix}.tbl.gz"), emit: tbl
    tuple val(meta), path("${prefix}/${prefix}.err.gz"), emit: err
    tuple val(meta), path("${prefix}/${prefix}.log.gz"), emit: log
    tuple val(meta), path("${prefix}/${prefix}.txt.gz"), emit: txt
    tuple val(meta), path("${prefix}/${prefix}.tsv.gz"), emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args   ?: ''
    prefix               = task.ext.prefix ?: "${meta.id}"
    def input            = fasta.toString() - ~/\.gz$/
    def decompress       = fasta.getExtension() == "gz" ? "gunzip -c ${fasta} > ${input}" : ""
    def cleanup          = fasta.getExtension() == "gz" ? "rm ${input}" : ""
    def proteins_opt     = proteins ? "--proteins ${proteins}" : ""
    def prodigal_tf_in   = prodigal_tf ? "--prodigaltf ${prodigal_tf}" : ""
    """
    ${decompress}

    prokka \\
        ${args} \\
        --cpus ${task.cpus} \\
        --prefix ${prefix} \\
        ${proteins_opt} \\
        ${prodigal_tf_in} \\
        ${input}

    ${cleanup}

    gzip ${prefix}/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.gff
    touch ${prefix}/${prefix}.gbk
    touch ${prefix}/${prefix}.fna
    touch ${prefix}/${prefix}.faa
    touch ${prefix}/${prefix}.ffn
    touch ${prefix}/${prefix}.sqn
    touch ${prefix}/${prefix}.fsa
    touch ${prefix}/${prefix}.tbl
    touch ${prefix}/${prefix}.err
    touch ${prefix}/${prefix}.log
    touch ${prefix}/${prefix}.txt
    touch ${prefix}/${prefix}.tsv
    touch ${prefix}/${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
