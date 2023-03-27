process FETCH_NCBI {
    tag "$ncbi_accs"
    label 'process_long'

    conda "anaconda::wget=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.20.1' :
        'quay.io/biocontainers/wget:1.20.1' }"
        container ""

    input:
    path ncbi_accs

    output:
    path "*.fna.gz"    , emit: fnas          // Contigs
    path "failures.txt", emit: failures
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir genomes/
    # Fetch file indexes
    wget -O 00refseq.index  ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
    wget -O 10genbank.index ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
    echo "Contig files for the following accessions were not found at NCBI" > failures.txt
    for a in \$(cat $ncbi_accs); do
        d=\$(grep \$a *.index | cut -f 20 | head -n 1 | sed 's/https/ftp/')
        echo "--> \$a: \$d <--"
        \$(wget \${d}/*_genomic.fna.gz)
        rm -f *from_genomic*.fna.gz
        if [ \$(ls *_genomic.fna.gz | wc -l) ]; then
            echo "\$a"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$( wget --version | grep '^GNU' | sed 's/GNU Wget //' | sed 's/ .*//' )
    END_VERSIONS
    """
}
