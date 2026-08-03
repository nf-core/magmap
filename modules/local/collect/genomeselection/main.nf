// Safely quote a Groovy value as a single-quoted R string literal. Used everywhere
// below that a sample ID or genome accession (documented as accepting arbitrary text
// via --genomeinfo) gets embedded into generated R source -- without this, a value
// containing a quote could break the generated R syntax, or worse.
def rq(v) {
    return "'" + v.toString().replace('\\', '\\\\').replace("'", "\\'") + "'"
}

process COLLECT_GENOMESELECTION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    tuple val(meta), val(genome_sets) // genome_sets: [ [ id: <sample or 'all'>, accnos: [ accno, ... ] ], ... ]
    val local_accnos                  // [ accno, ... ] -- genomes originating from --genomeinfo that were kept in the run

    output:
    tuple val(meta), path("${prefix}.genome_selection_mqc.tsv"), emit: multiqc
    tuple val(meta), path("${prefix}.genome_selection.tsv.gz") , emit: full_table
    path "versions.yml"                                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def local_accnos_r = local_accnos ?
        "c(${local_accnos.collect { a -> rq(a) }.join(', ')})" :
        "character(0)"
    def genome_set_rows = genome_sets.collectMany { s ->
        s.accnos.collect { accno -> "${rq(s.id)}, ${rq(accno)}," }
    }
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(readr)
    library(tidyr)

    local_accnos <- ${local_accnos_r}

    genome_sets <- tribble(
        ~sample, ~accno,
        ${genome_set_rows.join('\n        ')}
    ) %>%
        mutate(source = if_else(accno %in% local_accnos, 'local', 'remote'))

    write_tsv(genome_sets, "${prefix}.genome_selection.tsv.gz")

    summary_tab <- genome_sets %>%
        count(sample, source) %>%
        pivot_wider(names_from = source, values_from = n, values_fill = 0)
    if ( ! 'local'  %in% names(summary_tab) ) summary_tab\$local  <- 0
    if ( ! 'remote' %in% names(summary_tab) ) summary_tab\$remote <- 0
    summary_tab <- summary_tab %>%
        transmute(sample, n_local = local, n_remote = remote) %>%
        arrange(sample)

    writeLines(
        c(
            "# id: 'genome_selection'",
            "# section_name: 'Genome selection'",
            "# description: 'Number of local (user-provided) and remote (NCBI-fetched) genomes selected per genome set -- one row per sample when --genomeset_mode sample is used, a single row otherwise. See summary_tables/${prefix}.genome_selection.tsv.gz for the full per-genome breakdown.'",
            "# plot_type: 'table'"
        ),
        "${prefix}.genome_selection_mqc.tsv"
    )
    write_tsv(summary_tab, "${prefix}.genome_selection_mqc.tsv", append = TRUE, col_names = TRUE)

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    tidyr: ", packageVersion('tidyr'))
        ),
        "versions.yml"
    )
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genome_selection_mqc.tsv
    echo | gzip -c > ${prefix}.genome_selection.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        readr: 2.0.0
        tidyr: 1.1.3
    END_VERSIONS
    """
}
