process COLLECT_UNASSIGNEDCOUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9d1a8065dcebe35c513db5eb4a5237795aa4bae4756086f3c4319a820a8a647c/data' :
        'community.wave.seqera.io/library/custom_collectstats:c1f477e6251c36bb' }"

    input:
    tuple val(meta), path(summaries)

    output:
    tuple val(meta), path("${prefix}.Unassigned_NoFeatures.counts.tsv.gz"), path("${prefix}.Unassigned_Ambiguity.counts.tsv.gz"), emit: counts
    path "versions.yml"                                                                                                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(purrr)
    library(stringr)

    statuses <- c('Unassigned_NoFeatures', 'Unassigned_Ambiguity')

    d <- purrr::map_dfr(
        Sys.glob('*.featureCounts.tsv.summary'),
        function(file) {
            raw   <- read_tsv(file, show_col_types = FALSE)
            sample <- str_remove(names(raw)[2], '\\\\.sorted\\\\.bam\$')
            raw %>%
                rename(status = 1, count = 2) %>%
                filter(status %in% statuses) %>%
                transmute(sample, status, count)
        }
    )

    # accno/orf/chr/start/end/strand/length/tpm have no meaning for these pseudo-features,
    # but CUSTOM_COLLECTSTATS reads every fcs file together and needs matching columns,
    # so pad them out as NA here.
    d <- d %>%
        mutate(
            accno = NA_character_, orf = NA_character_, chr = NA_character_,
            start = NA_integer_, end = NA_integer_, strand = NA_character_, length = NA_integer_,
            tpm = NA_real_
        )

    d %>%
        filter(status == 'Unassigned_NoFeatures') %>%
        select(accno, orf, chr, start, end, strand, length, sample, count, tpm) %>%
        write_tsv("${prefix}.Unassigned_NoFeatures.counts.tsv.gz")

    d %>%
        filter(status == 'Unassigned_Ambiguity') %>%
        select(accno, orf, chr, start, end, strand, length, sample, count, tpm) %>%
        write_tsv("${prefix}.Unassigned_Ambiguity.counts.tsv.gz")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    purrr: ", packageVersion('purrr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip -c > ${prefix}.Unassigned_NoFeatures.counts.tsv.gz
    echo | gzip -c > ${prefix}.Unassigned_Ambiguity.counts.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        purrr: 0.3.4
        readr: 2.0.0
        stringr: 1.4.0
    END_VERSIONS
    """
}
