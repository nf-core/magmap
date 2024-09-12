//
// Create a BBMap index out of a set of fasta nucleotide files.
//
include { LOCAL_CAT as CAT_FNA  } from '../../modules/local/cat'
include { BBMAP_INDEX           } from '../../modules/nf-core/bbmap/index/main'

workflow CREATE_BBMAP_INDEX {
    take:
        ch_genome_fnas

    main:
        ch_versions = Channel.empty()
        CAT_FNA( [id:'cat.fna'], ch_genome_fnas.collect() )
        ch_versions = ch_versions.mix(CAT_FNA.out.versions)

        BBMAP_INDEX (CAT_FNA.out.concatenated_files.map{ it[1]})
        ch_versions = ch_versions.mix(BBMAP_INDEX.out.versions)

    emit:
    index         = BBMAP_INDEX.out.index.collect()
    genomes_fnas  = CAT_FNA.out.concatenated_files.collect()
    versions      = ch_versions
}
