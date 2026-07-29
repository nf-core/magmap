process DUCKDB_TSV2PARQUET {
    tag "summary tables"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/duckdb:1.5.5--b362367d86dce804'
        : 'community.wave.seqera.io/library/duckdb:1.5.5--09cfd9fc55b2e6d4'}"

    input:
    path tsvs

    output:
    path "*.parquet"    , emit: parquet
    path "versions.yml" , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import duckdb

    for f in "${tsvs}".split():
        out = f.removesuffix('.gz').removesuffix('.tsv') + '.parquet'
        duckdb.sql(
            f"COPY (SELECT * FROM read_csv('{f}', delim='\\t', header=true)) TO '{out}' (FORMAT PARQUET)"
        )

    with open('versions.yml', 'w') as fh:
        fh.write('"${task.process}":\\n')
        fh.write(f'    duckdb: {duckdb.__version__}\\n')
    """

    stub:
    """
    for f in ${tsvs}; do
        out=\$(basename "\$f" | sed 's/\\.tsv\\.gz\$//; s/\\.tsv\$//').parquet
        touch "\$out"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duckdb: 1.5.5
    END_VERSIONS
    """
}
