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
    path gtdb_metadata      // GTDB metadata files, species and completeness/contamination for remote genomes
    path checkm_metadata    // CheckM/CheckM2 output files, completeness/contamination for local genomes
    val  species_preference // 'local', 'completeness' or 'gtdb'

    output:
    path "${outfile_remote}", emit: kept_remote
    path "${outfile_local}" , emit: kept_local
    path "versions.yml"     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "magmap"
    outfile_remote = "${prefix}.kept_remote_accessions.txt"
    outfile_local  = "${prefix}.kept_local_accessions.txt"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(stringr)
    library(purrr)

    species_preference <- '${species_preference}'

    local_accessions  <- read_lines('${local_accessions}')
    remote_candidates <- read_lines('${remote_candidates}')

    local_species <- read_tsv(
        c('${gtdbtk_metadata.join('\', \'')}'),
        show_col_types = FALSE
    ) %>%
        transmute(accno = user_genome, species = str_extract(classification, '(?<=s__).+\$')) %>%
        filter(accno %in% local_accessions)

    remote_meta <- read_tsv(
        c('${gtdb_metadata.join('\', \'')}'),
        col_types = cols(
            .default = col_character(),
            checkm_completeness = col_double(), checkm_contamination = col_double()
        )
    ) %>%
        transmute(
            accno = str_remove(accession, '^.._'),
            species = str_extract(gtdb_taxonomy, '(?<=s__).+\$'),
            checkm_completeness, checkm_contamination
        ) %>%
        filter(accno %in% remote_candidates)

    if ( species_preference == 'local' ) {
        local_species_set <- local_species %>%
            filter(!is.na(species), species != '') %>%
            pull(species) %>%
            unique()

        # Drop a remote candidate only if it shares a classified species with a selected local
        # genome. Candidates missing from the GTDB metadata are kept (fail open): we can't
        # determine their species, so we don't drop a genome Sourmash asked for.
        dropped_remote <- remote_meta %>%
            filter(!is.na(species), species != '', species %in% local_species_set) %>%
            pull(accno)

        kept_remote <- setdiff(remote_candidates, dropped_remote)
        kept_local  <- local_accessions
    } else {
        # species_preference is 'completeness' or 'gtdb': pick the best-scoring genome per
        # species across BOTH local and remote candidates -- a local genome can lose here.
        checkm_cols <- c(
             '# contigs' = NA_integer_,
             'Genome size (bp)' = NA_integer_, 'Genome_Size' = NA_integer_,
             'Strain heterogeneity' = NA_real_,
             'contig_count' = NA_integer_
        )

        rcheckm <- function(fname) {
            read_tsv(fname, col_types = cols(.default = col_character())) %>%
                tibble::add_column(!!!checkm_cols[setdiff(names(checkm_cols), names(.))]) %>%
                    rename(accno = 1) %>%
                    transmute(
                        accno,
                        checkm_completeness = as.double(Completeness),
                        checkm_contamination = as.double(Contamination)
                    )
        }
        local_checkm <- purrr::map(c('${checkm_metadata.join('\', \'')}'), rcheckm) %>% purrr::reduce(union)

        candidates <- bind_rows(
            local_species %>% left_join(local_checkm, by = 'accno') %>% mutate(source = 'local'),
            remote_meta %>% mutate(source = 'remote')
        ) %>%
            mutate(
                score = if ( species_preference == 'gtdb' )
                    checkm_completeness - 5 * checkm_contamination
                else
                    checkm_completeness
            )

        # Fail open: keep anything we can't meaningfully compare (unclassified species or no score)
        unevaluable <- candidates %>% filter(is.na(species) | species == '' | is.na(score)) %>% pull(accno)

        best <- candidates %>%
            filter(!is.na(species), species != '', !is.na(score)) %>%
            group_by(species) %>%
            filter(score == max(score)) %>%
            # Ties go to the local genome, if one is among the tied candidates
            filter(!('local' %in% source) | source == 'local') %>%
            ungroup() %>%
            pull(accno)

        kept_accnos <- union(unevaluable, best)

        kept_remote <- intersect(remote_candidates, kept_accnos)
        kept_local  <- intersect(local_accessions, kept_accnos)
    }

    writeLines(kept_remote, "${outfile_remote}")
    writeLines(kept_local, "${outfile_local}")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    stringr: ", packageVersion('stringr')),
            paste0("    purrr: ", packageVersion('purrr'))
        ),
        "versions.yml"
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "magmap"
    outfile_remote = "${prefix}.kept_remote_accessions.txt"
    outfile_local  = "${prefix}.kept_local_accessions.txt"
    """
    echo $args

    cp ${remote_candidates} ${outfile_remote}
    cp ${local_accessions} ${outfile_local}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: 2.0.0
    END_VERSIONS
    """
}
