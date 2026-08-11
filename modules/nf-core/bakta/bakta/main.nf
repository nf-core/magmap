process BAKTA_BAKTA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50b75335f6394ae83fd05f364db27ee2eb75f4170e3525bb2aea47ad717a9e64/data'
        : 'community.wave.seqera.io/library/bakta_diamond:7830b94718da4f96'}"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf
    path regions
    path hmms

    output:
    tuple val(meta), path("${prefix}/${prefix}.embl.gz"), emit: embl
    tuple val(meta), path("${prefix}/${prefix}.faa.gz"), emit: faa
    tuple val(meta), path("${prefix}/${prefix}.ffn.gz"), emit: ffn
    tuple val(meta), path("${prefix}/${prefix}.fna.gz"), emit: fna
    tuple val(meta), path("${prefix}/${prefix}.gbff.gz"), emit: gbff
    tuple val(meta), path("${prefix}/${prefix}.gff3.gz"), emit: gff
    tuple val(meta), path("${prefix}/${prefix}.hypotheticals.tsv.gz"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}/${prefix}.hypotheticals.faa.gz"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}/${prefix}.tsv.gz"), emit: tsv
    tuple val(meta), path("${prefix}/${prefix}.txt.gz"), emit: txt
    tuple val(meta), path("${prefix}/${prefix}.json.gz"), emit: json

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf_opt = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    def regions_opt = regions ? "--regions ${regions}" : ""
    def hmms_opt = hmms ? "--hmms ${hmms}" : ""

    """
    ## Fake home due to fontconfig 'no writeable cache directory' issue
    mkdir nxf_home
    export HOME=\$PWD/nxf_home

    bakta \\
        ${fasta} \\
        ${args} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        --output ${prefix} \\
        ${proteins_opt} \\
        ${prodigal_tf_opt} \\
        ${regions_opt} \\
        ${hmms_opt} \\
        --db ${db}

    gzip ${prefix}/*
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export FONTCONFIG_PATH=\$PWD/.fontconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir .fontconfig .cache
    mkdir nxf_home
    export HOME=\$PWD/nxf_home

    mkdir ${prefix}
    touch ${prefix}/${prefix}.embl
    touch ${prefix}/${prefix}.faa
    touch ${prefix}/${prefix}.ffn
    touch ${prefix}/${prefix}.fna
    touch ${prefix}/${prefix}.gbff
    touch ${prefix}/${prefix}.gff3
    touch ${prefix}/${prefix}.hypotheticals.tsv
    touch ${prefix}/${prefix}.hypotheticals.faa
    touch ${prefix}/${prefix}.tsv
    touch ${prefix}/${prefix}.txt
    touch ${prefix}/${prefix}.json

    gzip ${prefix}/*
    """
}
