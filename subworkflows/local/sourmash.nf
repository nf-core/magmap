//
// select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER                   } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOMES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SAMPLES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'
include { FILTER_ACCNO                      } from '../../modules/local/create_accno_list'
include { COLLECTGENOMES                    } from '../../modules/local/collectgenomes'

workflow SOURMASH {
    take:
        reference_genomes
        samples_reads
        indexes
        reference_csv

    main:
        // I like that you create named variables for these, but they look more like config file
        // params than module arguments. Now, they *are* module arguments so we need to handle them,
        // but wouldn't it be better to let the subworkflow take them, and expose them via parameters?
        save_unassigned    = true
        save_matches_sig   = true
        save_prefetch      = true
        save_prefetch_csv  = false

        ch_versions = Channel.empty()

        GENOMES_SKETCH(reference_genomes)
        ch_versions = ch_versions.mix(GENOMES_SKETCH.out.versions)

        SAMPLES_SKETCH(samples_reads)
        ch_versions = ch_versions.mix(SAMPLES_SKETCH.out.versions)

        SAMPLES_SKETCH.out.signatures
            .collect{ it[1] }
            .map { [ [id: 'samples_sig'], it ] }
            .set { ch_sample_sigs }

        GENOMES_SKETCH.out.signatures
            .collect()
            .map { it[1] }
            .combine(indexes)
            .set { ch_genome_sigs }

        SOURMASH_GATHER ( ch_sample_sigs, ch_genome_sigs, save_unassigned, save_matches_sig, save_prefetch, save_prefetch_csv )
        ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

        SOURMASH_GATHER.out.result
            .map{ it[1] }
            .splitCsv( sep: ',', header: true, quote: '"')
            .map { it.name.replaceFirst(' .*', '') }
            .set { ch_accnos }

    COLLECTGENOMES(ch_accnos, reference_csv )

    emit:
        gindex        = GENOMES_SKETCH.out.signatures
        sindex        = SAMPLES_SKETCH.out.signatures
        accnos        = ch_accnos
        versions      = ch_versions
}
