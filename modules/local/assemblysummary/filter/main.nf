process ASSEMBLYSUMMARY_FILTER {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/0600fdb9287986d3d684b95d6e94f597c41b24a746608532d7c0836bd7b44c15/data' :
        'community.wave.seqera.io/library/gawk_gzip:fe0c202000bc8ed4' }"

    input:
    path accessions
    path summaries, stageAs: 'summaries/*'

    output:
    path "ncbi_genomeinfo.tsv", emit: tsv
    tuple val("${task.process}"), val('gawk'), eval("gawk 'BEGIN { print PROCINFO[\"version\"] }'"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf 'accno\\tftp_path\\n' > ncbi_genomeinfo.tsv
    gzip -cdf summaries/* | gawk -F'\\t' -v OFS='\\t' '
        NR == FNR { want[\$1]; next }
        /^#assembly_accession/ { for (i = 1; i <= NF; i++) if (\$i == "ftp_path") col = i; next }
        /^#/ { next }
        \$1 in want { print \$1, \$col }
    ' ${accessions} - >> ncbi_genomeinfo.tsv
    """

    stub:
    """
    printf 'accno\\tftp_path\\n' > ncbi_genomeinfo.tsv
    """
}
