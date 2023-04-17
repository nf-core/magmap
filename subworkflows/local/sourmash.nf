//
// select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER                   } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOMES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SAMPLES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'
include { GUNZIP                            } from '../../modules/nf-core/gunzip/main'
include { FILTER_ACCNO                      } from '../../modules/local/create_accno_list'

workflow SOURMASH {
    take:
        reference_genomes
        samples_reads
        indexes

    main:
        save_unassigned    = true
        save_matches_sig   = true
        save_prefetch      = true
        save_prefetch_csv  = false

        ch_versions = Channel.empty()

        GENOMES_SKETCH(reference_genomes)
        SAMPLES_SKETCH(samples_reads)
        ch_versions = ch_versions.mix(SAMPLES_SKETCH.out.versions)

        SOURMASH_GATHER(SAMPLES_SKETCH.out.signatures
                        .collect{ it[1] }
                        .map { [ [id: 'samples_sig'], it ] },
                        GENOMES_SKETCH.out.signatures
                        .collect()
                        .map { it[1] }
                        .combine(indexes),
                        save_unassigned,
                        save_matches_sig,
                        save_prefetch,
                        save_prefetch_csv
                        )

        GUNZIP(SOURMASH_GATHER.out.result.map{ [ it[0], it[1] ] } )

        FILTER_ACCNO( GUNZIP.out.gunzip)

    emit:
        gindex        = GENOMES_SKETCH.out.signatures
        sindex        = SAMPLES_SKETCH.out.signatures
        result        = FILTER_ACCNO.out.filtered_accno

        versions      = ch_versions

}
