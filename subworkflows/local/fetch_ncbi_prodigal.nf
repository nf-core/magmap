include { FETCH_NCBI } from '../../modules/local/fetch_ncbi'
include { GUNZIP     } from '../../modules/nf-core/gunzip/main'
include { PRODIGAL   } from '../../modules/nf-core/prodigal/main'

process TEST {
    input:
    tuple val(meta), path(file)
    val type

    output:
    path "${file}.out", emit: fsize

    script:
    """
    echo "$meta.id, $type" > ${file}.out
    ls -lL $file >> ${file}.out
    """
}

workflow FETCH_NCBI_PRODIGAL {
    take: ch_ncbi_accessions

    main:
    FETCH_NCBI(ch_ncbi_accessions)

    GUNZIP(FETCH_NCBI.out.fnas.flatten())

    //PRODIGAL(GUNZIP.out.gunzip.map { [ [ id: it.getBaseName() ], it ] }, 'gff')
    //TEST(GUNZIP.out.gunzip.map { [ id: 'hej' ], it })
    PRODIGAL(GUNZIP.out.gunzip.map { [ [ id: it.getBaseName() ], it ] }, 'gff')

    emit:
    fnas = FETCH_NCBI.out.fnas
    gffs = PRODIGAL.out.gene_annotations
}
