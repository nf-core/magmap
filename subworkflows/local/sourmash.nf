//
// select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER                   } '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOMES_SKETCH } '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SAMPLES_SKETCH } '../../modules/nf-core/sourmash/sketch/main'

workflow SOURMASH {
    take:
        reference_genomes
        samples_reads
    main:
        GENOMES_SKETCH(reference_genomes)
        SAMPLES_SKETCH(samples_reads)


}
