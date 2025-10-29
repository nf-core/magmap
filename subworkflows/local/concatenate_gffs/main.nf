//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CAT_MANY as CAT_GFF   } from '../../../modules/local/cat/many/main'
include { GENOMES2ORFS          } from '../../../modules/local/genomes_2_orfs'

workflow CONCATENATE_GFFS {
    take:
    ch_genome_gffs

    main:
        ch_versions = Channel.empty()

        CAT_GFF([id:'genomes'], ch_genome_gffs.collect())
        ch_versions = ch_versions.mix(CAT_GFF.out.versions)

        GENOMES2ORFS(ch_genome_gffs.collect().map { gffs -> [ [ id: 'genomes' ], gffs ] })
        ch_versions = ch_versions.mix(GENOMES2ORFS.out.versions)

    emit:
    gff          = CAT_GFF.out.concatenated_files
    genomes2orfs = GENOMES2ORFS.out.genomes2orfs
    versions     = ch_versions
}
