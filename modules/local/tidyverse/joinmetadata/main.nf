process TIDYVERSE_JOINMETADATA {
    tag "selected genomes"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546'
    }"

    input:
    path genomes            // Genome accessions to output information for
    path gtdb_metadata      // GTDB metadata files
    path gtdbtk_metadata    // Output files from GTDB-Tk
    path checkm_metadata    // Output files from CheckM/CheckM2

    output:
    path "*.genome_metadata.tsv.gz", emit: genome_metadata
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "magmap"
    def outfile = "${prefix}.genome_metadata.tsv.gz"

    def read_gtdb_metadata = """
        gtdb_metadata <- tibble(
            accno = character(), checkm_completeness = double(),
            checkm_contamination = double(), checkm_strain_heterogeneity = double(),
            contig_count = integer(), genome_size = integer(),
            gtdb_genome_representative = character(), gtdb_representative = logical(),
            gtdb_taxonomy = character(), source = character(), priority = integer()
        )
    """
    if ( gtdb_metadata ) {
        read_gtdb_metadata = """
            gtdb_metadata <- read_tsv(
                c('${gtdb_metadata.join('\', \'')}'),
                col_types = cols(
                    .default = col_character(), checkm_completeness = col_double(),
                    checkm_contamination = col_double(), checkm_strain_heterogeneity = col_double(),
                    contig_count = col_integer(), genome_size = col_integer(),
                    gtdb_representative = col_logical()
                )
            ) %>%
                transmute(
                    accno = str_remove(accession, '^.._'),
                    checkm_completeness, checkm_contamination, checkm_strain_heterogeneity,
                    contig_count, genome_size, gtdb_genome_representative, gtdb_representative,
                    gtdb_taxonomy, source = 'GTDB metadata', priority = 10
                )
        """
    }

    def read_gtdbtk_metadata = """
        gtdbtk_metadata <- tibble(accno = character(), gtdb_taxonomy = character())
    """
    if ( gtdbtk_metadata ) {
        read_gtdbtk_metadata = """
            gtdbtk_metadata <- read_tsv(c('${gtdbtk_metadata.join('\', \'')}'), show_col_types = FALSE) %>%
                transmute(accno = user_genome, gtdb_taxonomy = classification)
        """
    }

    def read_checkm_metadata = """
        checkm_metadata <- tibble(
            accno = character(), checkm_completeness = double(), checkm_contamination = double(),
            checkm_strain_heterogeneity = double(), contig_count = integer(),
            genome_size = integer()
        )
    """
    if ( checkm_metadata ) {
        read_checkm_metadata = """
        checkm_cols = c('# contigs' = NA_integer_, 'Genome size (bp)' = NA_integer_, 'Strain heterogeneity' = NA_integer_)
        checkm_metadata <- read_tsv(
            c('${checkm_metadata.join('\', \'')}'),
            col_types = cols(.default = col_character())
        ) %>%
            tibble::add_column(!!!checkm_cols[setdiff(names(checkm_cols), names(.))]) %>%
            rename(accno = 1) %>%
            transmute(
                accno,
                checkm_completeness = as.integer(Completeness), checkm_contamination = as.integer(Contamination),
                checkm_strain_heterogeneity = as.integer(`Strain heterogeneity`),
                contig_count = as.integer(`# contigs`), genome_size = as.integer(`Genome size (bp)`)
            )
        """
    }

    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    genomes <- read_tsv('${genomes}', col_types = 'c', col_names = c('accno'))

    ${read_gtdb_metadata}

    ${read_gtdbtk_metadata}

    ${read_checkm_metadata}

    genomes %>%
        left_join(
            gtdb_metadata %>%
                union(
                    gtdbtk_metadata %>% full_join(checkm_metadata, by = join_by(accno)) %>%
                        mutate(
                            gtdb_genome_representative = '', gtdb_representative = NA,
                            source = 'GTDB-Tk/CheckM', priority = 0
                        )
                ) %>%
                group_by(accno) %>%
                filter(priority == max(priority)) %>%
                ungroup(),
            by = join_by(accno)
        ) %>%
        arrange(accno) %>%
        write_tsv("${outfile}")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    tidyr: ", packageVersion('tidyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "magmap"
    """
    echo $args

    touch ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: 2.0.0
    END_VERSIONS
    """
}
