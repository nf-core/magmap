process GFFREAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff)
    path fasta

    output:
    tuple val(meta), path("*.gtf")  , emit: gtf             , optional: true
    tuple val(meta), path("*.gff3") , emit: gffread_gff     , optional: true
    tuple val(meta), path("*.fasta"), emit: gffread_fasta   , optional: true
    tuple val(meta), path("*.bed")  , emit: bed             , optional: true
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("--bed")    ? 'bed' : ( args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { flag -> args.contains(flag) } ) ? 'fasta' : 'gff3' ) )
    def gff_gz      = gff.name.endsWith('.gz')
    def fasta_gz    = fasta && fasta.name.endsWith('.gz')
    def fasta_in    = fasta_gz                  ? fasta.baseName : fasta
    def fasta_arg   = fasta                     ? "-g $fasta_in" : ''
    def output_name = "${prefix}.${extension}"
    def output      = extension == "fasta"      ? "$output_name" : "-o $output_name"
    def args_sorted = args.replaceAll(/(.*)(-[wxy])(.*)/) { _all, pre, param, post -> "$pre $post $param" }.trim()
    // args_sorted  = Move '-w', '-x', and '-y' to the end of the args string as gffread expects the file name after these parameters
    if ( "$output_name" in [ "$gff", "$fasta", "$fasta_in" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    // gffread cannot read gzipped input: a gzipped gff yields no features, without an error
    """
    ${fasta_gz ? "gunzip -c $fasta > $fasta_in" : ''}

    ${gff_gz ? "gunzip -c $gff |" : ''} gffread \\
        ${gff_gz ? '-' : gff} \\
        $fasta_arg \\
        $args_sorted \\
        $output

    # gffread names proteins by transcript (CDS parent); FeatureCounts uses the CDS ID
    ${gff_gz ? "gunzip -c $gff" : "cat $gff"} | awk -F '\\t' '
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
    ' - $output_name > renamed.tmp
    mv renamed.tmp $output_name

    ${fasta_gz ? "rm -f $fasta_in ${fasta_in}.fai" : ''}
    """

    stub:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("--bed")    ? 'bed' : ( args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { flag -> args.contains(flag) } ) ? 'fasta' : 'gff3' ) )
    def fasta_in    = fasta && fasta.name.endsWith('.gz') ? fasta.baseName : fasta
    def output_name = "${prefix}.${extension}"
    if ( "$output_name" in [ "$gff", "$fasta", "$fasta_in" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch $output_name
    """
}
