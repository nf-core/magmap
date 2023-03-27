process CONCATENATE {
    label 'process_long'

    conda (params.enable_conda ? "conda-forge::pigz=2.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0"
    }

    input:
    val  outfile
    path files

    output:
    path "${outfilename}", emit: file

    script:


    // If input is gzipped and options.args is set to something, i.e. filtering will be done, unzip input
    catcmd  = ( files[0] =~ /.gz$/ ) ? "unpigz -c -p" : "cat"

    // If input is gzipped, but not unzipped by catcmd, don't zip again
    outcmd  = ( files[0] =~ /.gz$/ ) ? '' : '| pigz -c -p'

    // Create a temporary file name to avoid collisions in a second call
    outfilename = ( outfile != '' ) ? outfile : File.createTempFile('outfile', '.gz').getName()

    """
    ${catcmd} $files ${options.args} ${outcmd} > $outfilename
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
