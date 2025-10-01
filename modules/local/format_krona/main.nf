process FORMAT_KRONA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546'
    }"

    input:
    tuple val(meta), path(krona_file)

    output:
    tuple val(meta), path("*.tsv"), emit: format_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    # Read krona file (count + tab-separated taxonomy)
    krona_data <- read_tsv("${krona_file}",
        col_names = c("count", "taxonomy"),
        col_types = "dc",
        show_col_types = FALSE
    )

    # Calculate total reads for fraction calculation
    total_reads <- sum(krona_data\$count, na.rm = TRUE)

    # Format taxonomy
    formatted_data <- krona_data %>%
        filter(count > 0) %>%
        mutate(
            fraction = count / total_reads,
            # remove prefixes (k__, p__, etc.)
            taxonomy = str_remove_all(taxonomy, "[a-z]__") %>%
                str_replace_all("_", " ")
        ) %>%
        separate(
            taxonomy,
            into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
            sep = "\t",
            fill = "right",   # fill missing ranks with NA
            remove = TRUE
        ) %>%
        replace_na(list(
            superkingdom = "unclassified", phylum = "unclassified", class = "unclassified",
            order = "unclassified", family = "unclassified", genus = "unclassified", species = "unclassified"
        )) %>%
        select(fraction, superkingdom, phylum, class, order, family, genus, species)

    # Write output
    write_tsv(formatted_data, "${prefix}.tsv")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    tidyverse: ", packageVersion('tidyverse'))
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat > ${prefix}.tsv << 'EOF'
    fraction	superkingdom	phylum	class	order	family	genus	species
    0.06582267833109018	Eukaryota	Chordata	Mammalia	Artiodactyla	Suidae	Sus	Sus scrofa
    0.0259084791386272	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	Escherichia coli
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.1
        tidyverse: 2.0.0
    END_VERSIONS
    """
}
