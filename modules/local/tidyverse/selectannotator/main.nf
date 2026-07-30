process TIDYVERSE_SELECTANNOTATOR {
    tag "genomes without a GFF"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    path no_gff_accessions // Single column file (no header): genomes lacking a GFF that need annotation
    path genome_metadata   // TIDYVERSE_JOINMETADATA output, incl. a gtdb_taxonomy column
    val  annotator         // 'prokka', 'bakta_supported_only' or 'bakta_all'

    output:
    path "${prefix}.bakta_accessions.txt"        , emit: bakta_accessions
    path "${prefix}.prokka_accessions.txt"       , emit: prokka_accessions
    path "${prefix}.unclassified_accessions.txt" , emit: unclassified_accessions
    path "versions.yml"                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "magmap"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(stringr)

    annotator <- '${annotator}'

    no_gff <- read_lines('${no_gff_accessions}')

    domains <- read_tsv('${genome_metadata}', col_types = cols(.default = col_character())) %>%
        transmute(accno, domain = str_extract(gtdb_taxonomy, '(?<=d__)[^;]+')) %>%
        filter(accno %in% no_gff, !is.na(domain), domain != '')

    if ( annotator == 'prokka' ) {
        bakta_accessions        <- character()
        prokka_accessions       <- no_gff
        unclassified_accessions <- character()
    } else if ( annotator == 'bakta_all' ) {
        bakta_accessions        <- no_gff
        prokka_accessions       <- character()
        unclassified_accessions <- character()
    } else {
        # bakta_supported_only: Bacteria to Bakta, everything else (incl. genomes we
        # couldn't classify by domain) falls back to Prokka.
        bakta_accessions        <- domains %>% filter(domain == 'Bacteria') %>% pull(accno)
        unclassified_accessions <- setdiff(no_gff, domains %>% pull(accno))
        prokka_accessions       <- setdiff(no_gff, bakta_accessions)
    }

    writeLines(bakta_accessions, "${prefix}.bakta_accessions.txt")
    writeLines(prokka_accessions, "${prefix}.prokka_accessions.txt")
    writeLines(unclassified_accessions, "${prefix}.unclassified_accessions.txt")

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
    prefix = task.ext.prefix ?: "magmap"
    """
    echo $args

    touch ${prefix}.bakta_accessions.txt
    touch ${prefix}.prokka_accessions.txt
    touch ${prefix}.unclassified_accessions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: 2.0.0
    END_VERSIONS
    """
}
