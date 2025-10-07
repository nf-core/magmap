//
// Select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER                   } from '../../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOME_SKETCH  } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_INDEX  as GENOME_INDEX   } from '../../../modules/nf-core/sourmash/index/main'
include { SOURMASH_SKETCH as SAMPLE_SKETCH  } from '../../../modules/nf-core/sourmash/sketch/main'
include { WGET as WGET_GENOME               } from '../../../modules/nf-core/wget/main'

workflow SOURMASH {
    take:
        ch_sample_reads             // Fastq files with reads for each sample [ val(meta), [ path(reads) ] ]
        ch_indexes                  // List of Sourmash indexs [ path(index) ]
        ch_user_genomeinfo          // User provided genomes [ path(genome) ]
        ch_remote_genome_sources    // Paths to genome information in NCBI format, i.e. containing at least the assembly_accession and ftp_path fields: path(csvfile)
        ksize                       // K-mere size to use: val(odd_int)
        save_unassigned             // Boolean value passed to sourmash/gather
        save_matches_sig            // Boolean value passed to sourmash/gather
        save_prefetch               // Boolean value passed to sourmash/gather
        save_prefetch_csv           // Boolean value passed to sourmash/gather

    main:
        ch_versions = Channel.empty()

        ch_ncbi_genomeinfo = ch_remote_genome_sources
                .splitCsv(skip: 1, header: true, sep: '\t')
                .map { row ->
                    [
                        accno: row["#assembly_accession"],
                        genome_fna: "${row.ftp_path}/${row.ftp_path - ~/\/$/ - ~/.*\//}_genomic.fna.gz",
                        genome_gff: ""
                    ]
                }

        GENOME_SKETCH(ch_user_genomeinfo.map { [ [ id: it.accno ], it.genome_fna ] })
        ch_versions = ch_versions.mix(GENOME_SKETCH.out.versions)

        SAMPLE_SKETCH(ch_sample_reads)
        ch_versions = ch_versions.mix(SAMPLE_SKETCH.out.versions)

        ch_sample_sigs = SAMPLE_SKETCH.out.signatures

        ch_genome_sigs = GENOME_SKETCH.out.signatures
            .collect { meta, sig -> sig }
            .map { sigs -> [ [ id: 'local-genomes' ], sigs ] }

        GENOME_INDEX(ch_genome_sigs, ksize)
        ch_versions = ch_versions.mix(GENOME_INDEX.out.versions)

        def i = 0
        ch_database = GENOME_INDEX.out.signature_index
            //.map { meta, sig -> sig }
            .mix(
                ch_indexes.map { index ->
                    [ [ id: sprintf("remoteidx_%02d", i++) ], index ]
                }
            )

        // To make sure that all combinations of sample signatures and indexes are gathered below,
        // combine the two channels.
        // (In theory, this should not be required as the command supposedly can take multiple samples
        // and multiple indexes. This does not return the full set of hits however.)
        ch_gather = ch_sample_sigs
            .combine(ch_database)

        SOURMASH_GATHER(ch_gather.map { it -> [ it[0], it[1] ] }, ch_gather.map { it -> [ it[2], it[3] ] }, save_unassigned, save_matches_sig, save_prefetch, save_prefetch_csv)
        ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

        // The genomes that were selected by sourmash can either be local genomes provided by the user
        // with --genomeinfo, or genomes we need to fetch from NCBI

        // 1. Find the genomes that were selected
        ch_genomes = ch_user_genomeinfo
            .map { genome -> [ [ genome.accno ], genome ] }
            .join(
                SOURMASH_GATHER.out.result
                    .map { meta, csv -> csv }
                    .splitCsv( sep: ',', header: true, quote: '"')
                    // Strip everything except accession number from NCBI-like names
                    .map { genome ->
                        def matcher = ( genome.name =~ /^(GC[A-Z]_[0-9]+\.[0-9]+)/ )
                        return matcher ?
                            [ [ matcher[0][0] ], [ accno: matcher[0][0] ] ] :
                            [ [ genome.name ], [ accno: genome.name ] ]
                    }
                    .unique(),
                remainder: true
            )
            .branch { genome ->
                local: genome[1] && genome[2]
                    return genome[1]
                ncbi:  genome[2]
                    return genome[2]
            }

        // 2. Fetch NCBI genomes
        WGET_GENOME(
            ch_genomes.ncbi
                .map { genome -> [ [ genome.accno ] ] }
                .join(
                    ch_ncbi_genomeinfo
                        .map { genome -> [ [ genome.accno ], genome ] }
                )
                .map { genome -> [ [ id: genome[1].accno ], genome[1].genome_fna ] }
        )
        ch_versions = ch_versions.mix(WGET_GENOME.out.versions)

        // 3. Mix the local and the newly fetched NCBI genomes
        ch_filtered_genomes = ch_genomes.local
            .mix(
                WGET_GENOME.out.outfile
                    .map { genome -> [ accno: genome[0].id, genome_fna: genome[1] ] }
            )

    emit:
        gindex           = GENOME_SKETCH.out.signatures
        sindex           = SAMPLE_SKETCH.out.signatures
        filtered_genomes = ch_filtered_genomes
        versions         = ch_versions
}
