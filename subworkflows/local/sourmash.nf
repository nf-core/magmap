//
// select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER                   } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOMES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_INDEX  as GENOMES_INDEX  } from '../../modules/nf-core/sourmash/index/main'
include { SOURMASH_SKETCH as SAMPLES_SKETCH } from '../../modules/nf-core/sourmash/sketch/main'

workflow SOURMASH {
    take:
        ch_samples_reads
        ch_indexes
        ch_user_genomeinfo
        ch_ncbi_genomeinfo_files

    main:
        // I like that you create named variables for these, but they look more like config file
        // params than module arguments. Now, they *are* module arguments so we need to handle them,
        // but wouldn't it be better to let the subworkflow take them, and expose them via parameters?
        save_unassigned    = true
        save_matches_sig   = true
        save_prefetch      = true
        save_prefetch_csv  = false

        ch_versions = Channel.empty()

        ch_ncbi_genomeinfo_files
                .splitCsv(sep: '\t')
                .map { file(it[0]) }
                .splitCsv(skip: 1, header: true, sep: '\t')
                .map {
                    [
                        accno: it["#assembly_accession"],
                        genome_fna: "${it.ftp_path}/${it.ftp_path - ~/\/$/ - ~/.*\//}_genomic.fna.gz",
                        genome_gff: ""
                    ]
                }
                .set { ch_ncbi_genomeinfo }

        GENOMES_SKETCH(ch_user_genomeinfo.map { [ [ id: it.accno ], it.genome_fna ] })
        ch_versions = ch_versions.mix(GENOMES_SKETCH.out.versions)

        SAMPLES_SKETCH(ch_samples_reads)
        ch_versions = ch_versions.mix(SAMPLES_SKETCH.out.versions)

        SAMPLES_SKETCH.out.signatures
            .collect{ it[1] }
            .map { [ [id: 'samples_sig'], it ] }
            .set { ch_sample_sigs }

        GENOMES_SKETCH.out.signatures
            .collect{ meta, sig -> [ sig ] }
            .map{ sig -> [ [id: 'signatures'], sig ] }
            .set { ch_genome_sigs }

        GENOMES_INDEX(ch_genome_sigs)

        GENOMES_INDEX.out.signature_index
            .map{ meta, sig -> sig }
            .mix(ch_indexes)
            .collect()
            .set{ ch_database }

        SOURMASH_GATHER ( SAMPLES_SKETCH.out.signatures, ch_database, save_unassigned, save_matches_sig, save_prefetch, save_prefetch_csv )
        ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

        SOURMASH_GATHER.out.result
            .map{ meta, csv -> csv }
            .splitCsv( sep: ',', header: true, quote: '"')
            .map { it.name }
            .unique()
            .set { ch_accnos }

        // Subset the two genome info channels to only contain those that Sourmash identified
        // The user supplied channel takes precedence, so start with that
        ch_accnos
            .join(ch_user_genomeinfo.map { [ it.accno, [ it ] ]} )
            .map { it[1][0] }
            .set{ ch_matching_user_genomes }

        ch_accnos
            .join(ch_matching_user_genomes.map { [ it.accno, [ it ] ]}, remainder: true)
            .filter{ !it[1] }
            .map{ it[0] }
            .join( ch_ncbi_genomeinfo.map { [ it.accno, [ it ] ] } )
            .map { it[1][0] }
            .map { [ accno: it.accno, genome_fna: file(it.genome_fna), genome_gff: ''] }
            .mix(ch_matching_user_genomes)
            .set { ch_filtered_genomes }

        // Add entries from NCBI for the ones missing in ch_filtered_genomes

    emit:
        gindex           = GENOMES_SKETCH.out.signatures
        sindex           = SAMPLES_SKETCH.out.signatures
        filtered_genomes = ch_filtered_genomes
        versions         = ch_versions
}
