process TIDYVERSE_JOINFEATURECOUNTSACCNO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    tuple val(meta), path(counts)
    path g2ids

    output:
    tuple val(meta), path("${prefix}.counts.tsv.gz"), emit: counts
    path "versions.yml"                             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)

    g2ids <- read_tsv("${g2ids}", col_types = 'ccc')

    read_tsv("${counts}", show_col_types = FALSE) %>%
        inner_join(g2ids %>% select(-genome), by = join_by(orf)) %>%
        relocate(accno) %>%
        write_tsv("${prefix}.counts.tsv.gz")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    dplyr: ", packageVersion('dplyr'))
        ),
        "versions.yml"
    )
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv
    gzip ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        readr: 2.0.0
        dplyr: 1.0.7
    END_VERSIONS
    """
}
