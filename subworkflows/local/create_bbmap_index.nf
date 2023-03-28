//
// Create a BBMap index out of a set of fasta nucleotide files.
//

include { CAT_CAT as FIRST_CAT     } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as SECOND_CAT    } from '../../modules/nf-core/cat/cat/main'
include { BBMAP_INDEX              } from '../../modules/nf-core/bbmap/index/main'

workflow CREATE_BBMAP_INDEX {
    take: ch_genome_fnas

    main:
        FIRST_CAT  (ch_genome_fnas)
        CAT_CAT.out.file_out.view()
        BBMAP_INDEX        (CAT_CAT.out.file_out.map{ it[1]})

    emit:
    index         = BBMAP_INDEX.out.index
    versions      = BBMAP_INDEX.out.versions
}
