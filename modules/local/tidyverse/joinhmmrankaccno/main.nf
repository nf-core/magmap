process TIDYVERSE_JOINHMMRANKACCNO {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data'
        : 'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546'}"

    input:
    tuple val(meta), path(hmmrank)
    path g2ids

    output:
    tuple val(meta), path("${prefix}.hmmrank.tsv.gz"), emit: hmmrank
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)

    g2ids <- read_tsv("${g2ids}", col_types = 'ccc')

    # Text columns keep hmmrank's number formatting; hmmrank's accno is the ORF, magmap's is the genome
    read_tsv("${hmmrank}", col_types = cols(.default = col_character())) %>%
        rename(orf = accno) %>%
        left_join(g2ids, by = join_by(orf)) %>%
        relocate(accno, genome, orf) %>%
        write_tsv("${prefix}.hmmrank.tsv.gz")

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
    touch ${prefix}.hmmrank.tsv
    gzip ${prefix}.hmmrank.tsv

    Rscript -e '
        writeLines(
            c(
                "\\"${task.process}\\":",
                paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                paste0("    readr: ", packageVersion("readr")),
                paste0("    dplyr: ", packageVersion("dplyr"))
            ),
            "versions.yml"
        )
    '
    """
}
