process COLLECT_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    tuple val(meta), path(inputfiles)
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

    library(data.table)
    library(dtplyr)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    setDTthreads($task.cpus)

    g2ids <- read_tsv("${g2ids}", col_types = 'ccc')

    tibble(f = Sys.glob('*.featureCounts.tsv')) %>%
        mutate(
            d = purrr::map(
                f,
                function(file) {
                    fread(file, sep = '\\t', skip = 1) %>%
                        melt(measure.vars = c(ncol(.)), variable.name = 'sample', value.name = 'count') %>%
                        lazy_dt() %>%
                        filter(count > 0) %>%
                        mutate(
                            sample = str_remove(sample, '.sorted.bam'),
                            r = count/Length
                        ) %>%
                        rename(orf = Geneid, chr = Chr, start = Start, end = End, strand = Strand, length = Length) %>%
                        group_by(sample) %>%
                        mutate(tpm = r/sum(r) * 1e6) %>% ungroup() %>%
                        select(-r) %>%
                        as_tibble()
                }
            )
        ) %>%
        tidyr::unnest(d) %>%
        select(-f) %>%
        inner_join(g2ids %>% select(-genome), by = join_by(orf)) %>%
        relocate(accno) %>%
        write_tsv("${prefix}.counts.tsv.gz")

        writeLines(
            c(
                "\\"${task.process}\\":",
                paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                paste0("    dplyr: ", packageVersion('dplyr')),
                paste0("    dtplyr: ", packageVersion('dtplyr')),
                paste0("    data.table: ", packageVersion('data.table'))
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
        dplyr: 1.0.7
        dtplyr: 1.1.0
        data.table: 1.14.0
    END_VERSIONS
    """
}
