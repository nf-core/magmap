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
        ksize

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

        GENOMES_INDEX(ch_genome_sigs, ksize)
        ch_versions = ch_versions.mix(GENOMES_INDEX.out.versions)

        GENOMES_INDEX.out.signature_index
            .map{ meta, sig -> sig }
            .mix( ch_indexes )
            .collect()
            .set{ ch_database }

        SOURMASH_GATHER ( ch_sample_sigs, ch_database, save_unassigned, save_matches_sig, save_prefetch, save_prefetch_csv )
        ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

        SOURMASH_GATHER.out.result
            .map { meta, csv -> csv }
            .splitCsv( sep: ',', header: true, quote: '"')
            .map { it.name }
            .unique()
            .map { name ->
                def matcher = (name =~ /(GCA_[0-9]+\.[0-9]+|GCF_[0-9]+\.[0-9]+)/)
                if (matcher) {
                    return matcher[0][0] // Return the matched pattern
                }
            }
            .set { ch_accnos_ncbi }

        SOURMASH_GATHER.out.result
            .map{ meta, csv -> csv }
            .splitCsv( sep: ',', header: true, quote: '"')
            .map { row -> row.name }
            .unique()
            .filter { name ->
                !(name =~ /(GCA_[0-9]+\.[0-9]+|GCF_[0-9]+\.[0-9]+)/)
            }
            .set { ch_all_non_ncbi_user_accnos }


        // Subset the two genome info channels to only contain those that Sourmash identified
        // The user supplied channel takes precedence, so start with that
        ch_all_non_ncbi_user_accnos
            .join(ch_user_genomeinfo.map { [ it.accno, [ it ] ]} )
            .map { it[1][0] }
            .set { ch_matching_user_non_ncbi_genomes }

        ch_accnos_ncbi
            .join(ch_user_genomeinfo.map { [ it.accno, [ it ] ]} )
            .map { it[1][0] }
            .set { ch_matching_user_ncbi_genomes }

        ch_accnos_ncbi
            .map { accno -> [accno, null] } // Initialize the channel with accno and null
            .join(ch_matching_user_ncbi_genomes.map { [it.accno, [ it ] ] }, remainder: true) // Perform the join
            .flatMap { tuple ->

                def accno = tuple[0] // accno from the left channel
                def matches = tuple[2] // Should be a list or null

                if (matches == null || matches.isEmpty()) {
                    return [[accno, null]]
                } else {
                    return matches.collect { match -> [accno, match] }
                }
            }
            .filter { it[1] == null } // Keep only tuples with null data
            .map { it[0] } // Extract accno
            .join(ch_ncbi_genomeinfo.map { [it.accno, it ] }) // Join with genome info
            .flatMap { tuple ->

                    def accno = tuple[0] // accno from the left channel
                    def genomeInfos = tuple[1] // Should be a list or null
                    if (genomeInfos == null || genomeInfos.isEmpty()) {
                        return [] // Discard tuples with null or empty genomeInfos
                    }
        
                    // Extract values from the map
                    def genomeFna = genomeInfos.genome_fna ?: ''
                    def genomeGff = genomeInfos.genome_gff ?: ''
    
                // Return the formatted result
                return [[accno: accno, genome_fna: genomeFna, genome_gff: genomeGff]]
            }
            .flatten()
            .mix(ch_matching_user_ncbi_genomes) // Ensure proper mixing with other data
            .mix(ch_matching_user_non_ncbi_genomes)
            .set { ch_filtered_genomes }

    emit:
        gindex           = GENOMES_SKETCH.out.signatures
        sindex           = SAMPLES_SKETCH.out.signatures
        filtered_genomes = ch_filtered_genomes
        versions         = ch_versions
}
