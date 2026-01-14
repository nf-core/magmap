process PROKKAGFF2TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546'
    }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${outfile}"), emit: tsv
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    outfile = "${prefix}/${prefix}.prokka-annotations.tsv.gz"

    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)

    dir.create('${prefix}')

    fread(
        cmd = "zgrep -P '\\t' $gff",
        col.names = c('contig', 'gene_caller', 'feature', 'start', 'end', 'a', 'strand', 'b', 'c')
    ) %>%
        separate_rows(c, sep = ';') %>%
        separate(c, c('k', 'v'), sep = '=') %>%
        pivot_wider(names_from = k, values_from = v) %>%
        select(-a, -b) %>%
        rename_all(str_to_lower) %>%
            mutate(
        db_xref   = if ("db_xref" %in% names(.)) db_xref else if ("dbxref" %in% names(.)) dbxref else NA_character_,
        gene      = if ("gene" %in% names(.)) gene else NA_character_,
        ec_number = if ("ec_number" %in% names(.)) ec_number else NA_character_
        ) %>%
        mutate(
            ec_number = if_else(
                is.na(ec_number) & str_detect(db_xref, 'EC:'),
                str_replace(db_xref, '.*EC:([0-9.-]+).*', '\\\\1'),
                ec_number
            )
        ) %>%
        transmute(
            locus_tag = id, ftype = feature, length = abs(start - end) + 1, gene, ec_number,
            cog = ifelse(str_detect(db_xref, 'COG'), str_replace(db_xref, '.*COG:(COG[0-9]+).*', '\\\\1'), NA),
            product
        ) %>%
        as.data.table() %>%
        write_tsv("${outfile}")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    data.table: ", packageVersion("data.table")),
            paste0("    dtplyr: "    , packageVersion("dtplyr")),
            paste0("    dplyr: "     , packageVersion("dplyr")),
            paste0("    tidyr: "     , packageVersion("tidyr")),
            paste0("    readr: "     , packageVersion("readr"))
        ),
        "versions.yml"
    )
    """
}
