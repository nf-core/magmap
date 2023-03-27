include { FETCH_NCBI } from '../../modules/local/fetch_ncbi'
include { GUNZIP     } from '../../modules/nf-core/gunzip/main'
include { PRODIGAL   } from '../../modules/nf-core/prodigal/main'

workflow FETCH_NCBI_PRODIGAL {
    take: ch_ncbi_accessions

    main:
    FETCH_NCBI(ch_ncbi_accessions)

    FETCH_NCBI.out.fnas.flatten().map {[ [ id: it.getBaseName() ], it ]}.set{ch_test}
    GUNZIP(ch_test)

    PRODIGAL(GUNZIP.out.gunzip, 'gff')

    emit:
    fnas = FETCH_NCBI.out.fnas
    gffs = PRODIGAL.out.gene_annotations
}
