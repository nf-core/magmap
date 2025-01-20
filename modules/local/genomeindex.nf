process GENOMEINDEX {
    label 'process_long'

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' }"

    input:
    path gffs

    output:
    path "${outfilename}", emit: genomes2id
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    outfilename = File.createTempFile('outfile', '.gz').getName()
    cpus        = Math.floor(task.cpus/2).toInteger()

    """
    echo -e "accno\tgenome\tID" | gzip -c > ${outfilename}

    for f in ${gffs}; do
        fn=\$(basename \$f .gff.gz)
        fn=\$(basename \$fn .gff)
        ac=\$(echo \$fn | sed 's/\\(G.._[0-9.]\\+\\)_.*/\\1/')

        if [ -f "\$f" ] && [ "\${f##*.}" = "gz" ]; then
            zcat "\$f"
        else
            cat "\$f"
        fi | grep -o 'ID=[A-Z0-9_]\\+' | \
            sed "s/^/\$ac\\t\$fn\\t/; s/ID=//" >> ${params.outfile_name}
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
