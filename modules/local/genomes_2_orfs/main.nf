process GENOMES2ORFS {
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'biocontainers/ubuntu:24.04' }"

    input:
    tuple val(meta), path(gffs, stageAs: 'gffs/*')

    output:
    path "*.genomes2orfs.tsv.gz", emit: genomes2orfs
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix  = task.ext.prefix ?: "${meta.id}"
    outfile = "${prefix}.genomes2orfs.tsv.gz"

    """
    echo -e "accno\tgenome\torf" | gzip -c > ${outfile}

    for f in ${gffs}; do
        fn=\$(basename \$f .gff.gz)
        fn=\$(basename \$fn .gff)
        ac=\$(echo \$fn | sed 's/\\(G.._[0-9.]\\+\\)_.*/\\1/')

        if [ -f "\$f" ] && [ "\${f##*.}" = "gz" ]; then
            zcat "\$f"
        else
            cat "\$f"
        fi | grep -o 'ID=[A-Z0-9_]\\+' | \
            sed "s/^/\$ac\\t\$fn\\t/; s/ID=//" | gzip -c >> ${outfile}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | grep "^sed" | sed 's/^.* //')
    END_VERSIONS
    """
}
