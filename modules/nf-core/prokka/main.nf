process PROKKA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl5321hdfd78af_4' :
        'biocontainers/prokka:1.14.6--pl5321hdfd78af_4' }"

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}/*.gff.gz"), emit: gff
    tuple val(meta), path("${prefix}/*.gbk.gz"), emit: gbk
    tuple val(meta), path("${prefix}/*.fna.gz"), emit: fna
    tuple val(meta), path("${prefix}/*.faa.gz"), emit: faa
    tuple val(meta), path("${prefix}/*.ffn.gz"), emit: ffn
    tuple val(meta), path("${prefix}/*.sqn.gz"), emit: sqn
    tuple val(meta), path("${prefix}/*.fsa.gz"), emit: fsa
    tuple val(meta), path("${prefix}/*.tbl.gz"), emit: tbl
    tuple val(meta), path("${prefix}/*.err.gz"), emit: err
    tuple val(meta), path("${prefix}/*.log.gz"), emit: log
    tuple val(meta), path("${prefix}/*.txt.gz"), emit: txt
    tuple val(meta), path("${prefix}/*.tsv.gz"), emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    """
    prokka \\
        $args \\
        --cpus $task.cpus \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf \\
        $fasta
    
    gzip $prefix/* 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
