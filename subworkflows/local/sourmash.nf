//
// select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER                   } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOMES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SAMPLES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'

workflow SOURMASH {
    take:
        reference_genomes
        samples_reads

    main:
        save_unassigned    = true
        save_matches_sig   = true
        save_prefetch      = true
        save_prefetch_csv  = true

        ch_versions = Channel.empty()

        GENOMES_SKETCH(reference_genomes)
        SAMPLES_SKETCH(samples_reads)
        ch_versions = ch_versions.mix(SAMPLES_SKETCH.out.versions)

        SOURMASH_GATHER(SAMPLES_SKETCH.out.signatures,
                        GENOMES_SKETCH.out.signatures.collect()
                        .map { it[1] },
                        save_unassigned,
                        save_matches_sig,
                        save_prefetch,
                        save_prefetch_csv
                        )

    emit:
        gindex        = GENOMES_SKETCH.out.signatures
        sindex        = SAMPLES_SKETCH.out.signatures
        versions      = ch_versions

}
