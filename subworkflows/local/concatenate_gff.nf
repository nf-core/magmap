//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CAT_CAT as FIRST_CAT  } from '../../modules/nf-core/cat/cat'
include { CAT_CAT as SECOND_CAT } from '../../modules/nf-core/cat/cat'

workflow CAT_GFFS {
    take: ch_genome_gffs

    main:
        FIRST_CAT   ('', ch_genome_gffs.collate(1000))
        SECOND_CAT  ('genomes.gff.gz', FIRST_CONCATENATE.out.file.collect())

    emit:
    gff    = SECOND_CAT.out.file
}
