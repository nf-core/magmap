process CHECK_DUPLICATES {
    label 'process_low'

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'biocontainers/ubuntu:24.04' }"

    input:
    path fnas

    output:
    path "duplicate-contigs.txt" , emit: duplicate_list, optional: true
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    # Try to find and count duplicate contig headers
    zgrep -h '>' *.fna.gz | sort | uniq -c | grep -v ' 1 ' > ./duplicate-contigs.txt || true

    # Check if duplicate-contigs.txt has any real content
    if [ \$(wc -l < ./duplicate-contigs.txt) -gt 0 ]; then
        echo "You have duplicate contig names"
        cat ./duplicate-contigs.txt
        exit 2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zgrep: \$( zgrep --version | sed 's/.*/1.5/' | head -n 1 )
    END_VERSIONS
    """
}
