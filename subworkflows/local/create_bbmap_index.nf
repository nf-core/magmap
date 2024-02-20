//
// Create a BBMap index out of a set of fasta nucleotide files.
//

include { CAT_CAT as FIRST_CAT     } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as SECOND_CAT    } from '../../modules/nf-core/cat/cat/main'
include { BBMAP_INDEX              } from '../../modules/nf-core/bbmap/index/main'

workflow CREATE_BBMAP_INDEX {
    take:
        ch_genome_fnas

    main:
        ch_versions = Channel.empty()
        FIRST_CAT   (ch_genome_fnas)
        ch_versions = ch_versions.mix(FIRST_CAT.out.versions)

        SECOND_CAT  (FIRST_CAT.out.file_out
                        .collect{ it[1] }
                        .map{ [ [ id:'references_fnas' ], it ] } )
        ch_versions = ch_versions.mix(SECOND_CAT.out.versions)

        BBMAP_INDEX (SECOND_CAT.out.file_out.map{ it[1]})
        ch_versions = ch_versions.mix(BBMAP_INDEX.out.versions)

    emit:
    index         = BBMAP_INDEX.out.index
    genomes_fnas  = SECOND_CAT.out.file_out
    versions      = ch_versions
}
