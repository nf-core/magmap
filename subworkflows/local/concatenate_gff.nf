//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CAT_MANY as CAT_GFF   } from '../../modules/local/cat_many'
include { GENOMEINDEX           } from '../../modules/local/genomeindex'
include { CAT_CAT as GINDEX_CAT } from '../../modules/nf-core/cat/cat'
include { PROKKAGFF2TSV         } from '../../modules/local/prokkagff2tsv'

workflow CAT_GFFS {
    take: ch_genome_gffs

    main:
        ch_versions = Channel.empty()

        CAT_GFF( [id:'cat.gff'], ch_genome_gffs.collect() )
        ch_versions = ch_versions.mix(CAT_GFF.out.versions)

        GENOMEINDEX(ch_genome_gffs.collect())
        ch_versions = ch_versions.mix(GENOMEINDEX.out.versions)

        GINDEX_CAT(GENOMEINDEX.out.genomes2id.collect().map { [ [id: 'genomes_index'], it ] })
        ch_versions = ch_versions.mix(GINDEX_CAT.out.versions)

        PROKKAGFF2TSV(CAT_GFF.out.concatenated_files)
        ch_versions = ch_versions.mix(PROKKAGFF2TSV.out.versions)

    emit:
    gff      = CAT_GFF.out.concatenated_files
    gindex   = GINDEX_CAT.out.file_out
    gfftsv   = PROKKAGFF2TSV.out.tsv
    versions = ch_versions
}
