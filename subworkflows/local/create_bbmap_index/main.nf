//
// Create a BBMap index out of a set of fasta nucleotide files.
//
include { CAT_MANY as CAT_FNA } from '../../../modules/local/cat/many/main'
include { BBMAP_INDEX         } from '../../../modules/nf-core/bbmap/index/main'

workflow CREATE_BBMAP_INDEX {
    take:
        ch_genome_fnas

    main:
        CAT_FNA( [id:'cat.fna'], ch_genome_fnas.collect() )

        BBMAP_INDEX (CAT_FNA.out.concatenated_files.map{ it -> it[1]})

    emit:
    index         = BBMAP_INDEX.out.index.collect()
    genome_fnas   = CAT_FNA.out.concatenated_files.collect()
}
