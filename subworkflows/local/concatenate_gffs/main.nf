//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CAT_MANY as CAT_GFF   } from '../../../modules/local/cat/many/main'
include { GENOMES2ORFS          } from '../../../modules/local/genomes_2_orfs'

workflow CONCATENATE_GFFS {
    take:
    ch_genome_gffs

    main:
        CAT_GFF([id:'genomes'], ch_genome_gffs.collect())

        GENOMES2ORFS(ch_genome_gffs.collect().map { gffs -> [ [ id: 'genomes' ], gffs ] })

    emit:
    gff          = CAT_GFF.out.concatenated_files
    genomes2orfs = GENOMES2ORFS.out.genomes2orfs
}
