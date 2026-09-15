// Safely quote a Groovy value as a single-quoted R string literal -- meta.id only has to
// match /^\S+$/, so a sample name with a quote character would otherwise break R syntax.
def rq(v) {
    return "'" + v.toString().replace('\\', '\\\\').replace("'", "\\'") + "'"
}

process TIDYVERSE_SPLITFEATURECOUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9d1a8065dcebe35c513db5eb4a5237795aa4bae4756086f3c4319a820a8a647c/data' :
        'community.wave.seqera.io/library/custom_collectstats:c1f477e6251c36bb' }"

    input:
    tuple val(meta), path(counts)
    path ftypes

    output:
    tuple val(meta), path("${prefix}.*.featureCounts.tsv"), emit: counts
    path "versions.yml"                                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)

    ftypes <- read_tsv("${ftypes}", col_types = cols(orf = col_character(), .default = col_guess())) %>%
        select(orf, ftype)

    counts <- read_tsv("${counts}", skip = 1, col_types = cols(Geneid = col_character(), .default = col_guess()))

    # Reconstructs the per-feature-type files a dedicated per-type featureCounts run
    # would have produced, so CUSTOM_COLLECTFEATURECOUNTS/TIDYVERSE_JOINFEATURECOUNTSACCNO
    # stay unchanged even though FEATURECOUNTS now runs once per sample.
    annotated <- counts %>%
        inner_join(ftypes, by = c('Geneid' = 'orf'))

    # Iterate over the requested feature types (meta.feature), not unique(annotated\$ftype):
    # a type with zero matching rows must still produce a (header-only) file, since the
    # output below is a mandatory glob requiring at least one match.
    for ( f in strsplit(${rq(meta.feature)}, ',')[[1]] ) {
        outfile <- paste0(${rq(prefix)}, ".", f, ".featureCounts.tsv")
        writeLines("# Split by TIDYVERSE_SPLITFEATURECOUNTS", outfile)
        annotated %>%
            filter(ftype == f) %>%
            select(-ftype) %>%
            write_tsv(outfile, append = TRUE, col_names = TRUE)
    }

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr'))
        ),
        "versions.yml"
    )
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Single-quoted and quote-escaped: meta.id only has to match /^\S+$/, so a shell
    // metacharacter here would otherwise break out of the unquoted touch command.
    safePrefix = prefix.replace("'", "'\\''")
    """
    touch '${safePrefix}.CDS.featureCounts.tsv'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        readr: 2.0.0
    END_VERSIONS
    """
}
