process GFFREAD_PROTEINS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4':
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(fna), path(gff)

    output:
    tuple val(meta), path("${prefix}.faa.gz"), emit: faa
    tuple val("${task.process}"), val('gffread'), eval("gffread --version 2>&1"), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    ${fna.name.endsWith('.gz') ? 'gunzip -c' : 'cat'} ${fna} > genome.fna
    ${gff.name.endsWith('.gz') ? 'gunzip -c' : 'cat'} ${gff} > genome.gff

    gffread genome.gff -g genome.fna ${args} -y proteins.faa

    # gffread names proteins by transcript (CDS parent); FeatureCounts uses the CDS ID
    awk -F '\\t' '
        FNR == NR {
            if (\$3 != "CDS") next
            id = ""
            parent = ""
            n = split(\$9, attrs, ";")
            for (i = 1; i <= n; i++) {
                if (attrs[i] ~ /^ID=/) id = substr(attrs[i], 4)
                else if (attrs[i] ~ /^Parent=/) parent = substr(attrs[i], 8)
            }
            if (parent == "") parent = id
            if (!(parent in cds)) cds[parent] = id
            next
        }
        /^>/ {
            name = substr(\$1, 2)
            if (name in cds) \$0 = ">" cds[name]
        }
        { print }
    ' genome.gff proteins.faa | gzip -c > ${prefix}.faa.gz

    rm genome.fna genome.fna.fai genome.gff proteins.faa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.faa.gz
    """
}
