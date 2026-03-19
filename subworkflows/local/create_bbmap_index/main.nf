//
// Create a BBMap index out of a set of fasta nucleotide files.
//
include { CAT_MANY as CAT_FNA } from '../../../modules/local/cat/many/main'
include { BBMAP_INDEX         } from '../../../modules/nf-core/bbmap/index/main'

workflow CREATE_BBMAP_INDEX {
    take:
        ch_genome_fnas // [ meta, list_of_fnas ]

    main:
        CAT_FNA(ch_genome_fnas.map { g -> [ id: g[0].id ] }, ch_genome_fnas.map { g -> g[1] })

        BBMAP_INDEX(CAT_FNA.out.concatenated_files)

    emit:
    index         = BBMAP_INDEX.out.index
    genome_fnas   = CAT_FNA.out.concatenated_files
}
