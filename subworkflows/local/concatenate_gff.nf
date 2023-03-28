//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CAT_CAT as FIRST_CAT  } from '../../modules/nf-core/cat/cat'
include { CAT_CAT as SECOND_CAT } from '../../modules/nf-core/cat/cat'
include { GENOMEINDEX           } from '../../modules/local/genomeindex'
include { CAT_CAT as GINDEX_CAT } from '../../modules/nf-core/cat/cat'

workflow CAT_GFFS {
    take: ch_genome_gffs

    main:
        ch_versions = Channel.empty()

        ch_genome_gffs
            .collect{ it[2] }
            .collate(1000)
            .map{ [ [ id: 'all_references'], it[0] ] }
            .set { ch_reference_gffs }

        FIRST_CAT   (ch_reference_gffs)
        SECOND_CAT  (FIRST_CAT.out.file_out.map{ [ [ id:'reference_gffs' ], it[1] ] })
        ch_versions = ch_versions.mix(SECOND_CAT.out.versions)

        GENOMEINDEX(ch_genome_gffs.collect{ it[2] })
        GINDEX_CAT(GENOMEINDEX.out.genomes2id.collect().map { [ [id: 'genomes_index'], it ] })

    emit:
    gff      = SECOND_CAT.out.file_out
    gindex   = GINDEX_CAT.out.file_out
    versions = ch_versions
}
