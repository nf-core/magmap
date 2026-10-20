// Quote a value as a single-quoted SQL string literal.
def sqlq(v) {
    return "'" + v.toString().replace("'", "''") + "'"
}

// Quote a value as a single-quoted shell word; meta.id may contain quote characters.
def shq(v) {
    return "'" + v.toString().replace("'", "'\\''") + "'"
}

process DUCKDB_SPLITFEATURECOUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68261e24307fdf80b9988d6fd13cf3551735e6dc0e7e38003259f7d8efa84cdb/data' :
        'community.wave.seqera.io/library/duckdb-cli:1.5.5--9c6d18d9f687a45d' }"

    input:
    tuple val(meta), path(counts)
    path ftypes

    output:
    tuple val(meta), path("${prefix}.*.featureCounts.tsv"), emit: counts
    tuple val("${task.process}"), val('duckdb'), eval("duckdb -noheader -list -c 'SELECT substr(version(), 2)'"), topic: versions, emit: versions_duckdb

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // One query per feature type: the build side of each join holds only that type's ORFs.
    def splits = meta.feature.tokenize(',').collect { f ->
        def outfile = "${prefix}.${f}.featureCounts.tsv"
        def part    = "${prefix}.${f}.part"
        def sql = """\
            SET threads=${task.cpus};
            SET memory_limit='${(task.memory.toMega() * 0.8) as long}MB';
            SET temp_directory='.';
            COPY (
                SELECT c.* EXCLUDE (rn)
                FROM (
                    SELECT *, row_number() OVER () AS rn
                    FROM read_csv(${sqlq(counts)}, delim='\\t', header=true, skip=1, all_varchar=true)
                ) c
                JOIN (
                    SELECT orf FROM read_csv(${sqlq(ftypes)}, delim='\\t', header=true, all_varchar=true)
                    WHERE ftype = ${sqlq(f)}
                ) t ON c.Geneid = t.orf
                ORDER BY c.rn
            ) TO ${sqlq(part)} (FORMAT csv, DELIMITER '\\t', HEADER);
"""
        // The quoted heredoc body is not expanded by the shell; its delimiter is the body's own hash.
        def tag = "SQL_${sql.md5()}"
        """
        duckdb <<'${tag}'
${sql}${tag}
        printf '%s\\n' '# Split by DUCKDB_SPLITFEATURECOUNTS' | cat - ${shq(part)} > ${shq(outfile)}
        rm ${shq(part)}
        """
    }
    splits.join('\n')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${shq("${prefix}.CDS.featureCounts.tsv")}
    """
}
