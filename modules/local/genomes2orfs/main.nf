process GENOMES2ORFS {
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'biocontainers/ubuntu:24.04' }"

    input:
    tuple val(meta), path(gffs, stageAs: 'gffs/*')

    output:
    tuple val(meta), path(outfile), emit: genomes2orfs
    tuple val("${task.process}"), val('sed'), eval('sed --version 2>&1 | grep "^sed" | sed "s/^.* //"'), emit: versions_sed, topic: versions

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
        fi | grep -o 'ID=[A-Z_0-9]\\+' | \
            sed "s/^/\$ac\\t\$fn\\t/; s/ID=//" | gzip -c >> ${outfile}
    done
    """

    stub:
    prefix  = task.ext.prefix ?: "${meta.id}"
    outfile = "${prefix}.genomes2orfs.tsv.gz"
    """
    echo "" | gzip -c > ${outfile}
    """
}
