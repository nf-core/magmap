//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CAT_CAT as FIRST_CAT  } from '../../modules/nf-core/cat/cat'
include { CAT_CAT as SECOND_CAT } from '../../modules/nf-core/cat/cat'
include { GENOMEINDEX           } from '../../modules/local/genomeindex'
include { CAT_CAT as GINDEX_CAT } from '../../modules/nf-core/cat/cat'
include { PROKKAGFF2TSV         } from '../../modules/local/prokkagff2tsv'

workflow CAT_GFFS {
    take: ch_genome_gffs

    main:
        ch_versions = Channel.empty()
        def i = 0
        ch_genome_gffs
            .map{ it[1] }
            .flatten()
            .collate(1000)
            .map{ [ [ id: "all_references${i++}" ], it ] }
            .set { ch_reference_gffs }
        FIRST_CAT   (ch_reference_gffs)
        ch_versions = ch_versions.mix(FIRST_CAT.out.versions)

        SECOND_CAT  (FIRST_CAT.out.file_out.map{ [ [ id:'reference_gffs' ], it[1] ] })
        ch_versions = ch_versions.mix(SECOND_CAT.out.versions)

        GENOMEINDEX(ch_genome_gffs.collect{ it[1] })
        ch_versions = ch_versions.mix(GENOMEINDEX.out.versions)

        GINDEX_CAT(GENOMEINDEX.out.genomes2id.collect().map { [ [id: 'genomes_index'], it ] })
        ch_versions = ch_versions.mix(GINDEX_CAT.out.versions)

        PROKKAGFF2TSV(SECOND_CAT.out.file_out)
        ch_versions = ch_versions.mix(PROKKAGFF2TSV.out.versions)

    emit:
    gff      = SECOND_CAT.out.file_out
    gindex   = GINDEX_CAT.out.file_out
    gfftsv   = PROKKAGFF2TSV.out.tsv
    versions = ch_versions
}
