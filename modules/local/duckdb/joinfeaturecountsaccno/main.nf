// Quote a value as a single-quoted SQL string literal.
def sqlq(v) {
    return "'" + v.toString().replace("'", "''") + "'"
}

// Quote a value as a single-quoted shell word; meta.id may contain quote characters.
def shq(v) {
    return "'" + v.toString().replace("'", "'\\''") + "'"
}

process DUCKDB_JOINFEATURECOUNTSACCNO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68261e24307fdf80b9988d6fd13cf3551735e6dc0e7e38003259f7d8efa84cdb/data' :
        'community.wave.seqera.io/library/duckdb-cli:1.5.5--9c6d18d9f687a45d' }"

    input:
    tuple val(meta), path(counts)
    path g2ids

    output:
    tuple val(meta), path("${prefix}.counts.tsv.gz"), emit: counts
    tuple val("${task.process}"), val('duckdb'), eval("duckdb -noheader -list -c 'SELECT substr(version(), 2)'"), topic: versions, emit: versions_duckdb

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def sql = """\
        SET threads=${task.cpus};
        SET memory_limit='${(task.memory.toMega() * 0.8) as long}MB';
        SET temp_directory='.';
        COPY (
            SELECT g.accno, c.* EXCLUDE (rn)
            FROM (
                SELECT *, row_number() OVER () AS rn
                FROM read_csv(${sqlq(counts)}, delim='\\t', header=true, all_varchar=true)
            ) c
            JOIN (
                SELECT accno, orf FROM read_csv(${sqlq(g2ids)}, delim='\\t', header=true, all_varchar=true)
            ) g ON c.orf = g.orf
            ORDER BY c.rn
        ) TO ${sqlq("${prefix}.counts.tsv.gz")} (FORMAT csv, DELIMITER '\\t', HEADER);
"""
    // The quoted heredoc body is not expanded by the shell; its delimiter is the body's own hash.
    def tag = "SQL_${sql.md5()}"
    """
    duckdb <<'${tag}'
${sql}${tag}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${shq("${prefix}.counts.tsv")}
    gzip ${shq("${prefix}.counts.tsv")}
    """
}
