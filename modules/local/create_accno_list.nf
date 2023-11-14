process FILTER_ACCNO {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' }"

    input:
    tuple val(meta), path(sourmashoutput)

    output:
    tuple val(meta), path("*.tsv.gz")   , emit: filtered_accno
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(stringr)

    sourmash_output <- read_csv('samples_sig.csv')

    filtered_genomes <- sourmash_output[,10] %>%
        as_data_frame() %>%
        separate(name, c('new','col'), sep = ',') %>%
        separate(new, c('accno','x'), sep = ' ') %>%
        as_tibble() %>%
        select(-x,-col) %>%
        write_tsv("${prefix}_filtered_accno.tsv.gz")

        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")), paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    dtplyr: ", packageVersion('dtplyr')), paste0("    data.table: ", packageVersion('data.table')) ), "versions.yml")

        """
}
