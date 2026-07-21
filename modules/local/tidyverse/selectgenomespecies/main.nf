process TIDYVERSE_SELECTGENOMESPECIES {
    tag "remote genomes"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    path local_accessions   // One accno per line: local genomes selected for this run
    path remote_candidates  // One accno per line: remote genomes Sourmash wants to fetch
    path gtdbtk_metadata    // GTDB-Tk output files, species assignments for local genomes
    path gtdb_metadata      // GTDB metadata files, species assignments for remote genomes

    output:
    path "${outfile}"  , emit: kept
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "magmap"
    outfile = "${prefix}.kept_remote_accessions.txt"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(stringr)

    local_accessions  <- read_lines('${local_accessions}')
    remote_candidates <- read_lines('${remote_candidates}')

    local_species <- read_tsv(
        c('${gtdbtk_metadata.join('\', \'')}'),
        show_col_types = FALSE
    ) %>%
        transmute(accno = user_genome, species = str_extract(classification, '(?<=s__).+\$')) %>%
        filter(accno %in% local_accessions, !is.na(species), species != '') %>%
        pull(species) %>%
        unique()

    remote_species <- read_tsv(
        c('${gtdb_metadata.join('\', \'')}'),
        col_types = cols(.default = col_character())
    ) %>%
        transmute(accno = str_remove(accession, '^.._'), species = str_extract(gtdb_taxonomy, '(?<=s__).+\$')) %>%
        filter(accno %in% remote_candidates)

    # Drop a remote candidate only if it shares a classified species with a selected local
    # genome. Candidates missing from the GTDB metadata are kept (fail open): we can't
    # determine their species, so we don't drop a genome Sourmash asked for.
    dropped <- remote_species %>%
        filter(!is.na(species), species != '', species %in% local_species) %>%
        pull(accno)

    kept <- setdiff(remote_candidates, dropped)

    writeLines(kept, "${outfile}")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "magmap"
    outfile = "${prefix}.kept_remote_accessions.txt"
    """
    echo $args

    cp ${remote_candidates} ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: 2.0.0
    END_VERSIONS
    """
}
