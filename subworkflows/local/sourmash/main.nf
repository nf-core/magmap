//
// Select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER as GATHER_USER_GENOMES   } from '../../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_GATHER as GATHER_REMOTE_GENOMES } from '../../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOME_SKETCH         } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_INDEX  as GENOME_INDEX          } from '../../../modules/nf-core/sourmash/index/main'
include { SOURMASH_SKETCH as SAMPLE_SKETCH         } from '../../../modules/nf-core/sourmash/sketch/main'
include { WGET as WGET_GENOME                      } from '../../../modules/nf-core/wget/main'

workflow SOURMASH {
    take:
        ch_sample_reads             // Fastq files with reads for each sample [ val(meta), [ path(reads) ] ]
        ch_indexes                  // List of Sourmash indexs [ path(index) ]
        index_list                  // Value of the indexes param, used for if clauses
        ch_user_genomeinfo          // User provided genomes [ path(genome) ]
        ch_remote_genome_sources    // Paths to genome information in NCBI format, i.e. containing at least the assembly_accession and ftp_path fields: path(csvfile)
        ksize                       // K-mere size to use: val(odd_int)
        skip_sourmash               // Boolean that controls whether user-provided genomes are sketched, indexed and used in gathering genomes

    main:
        ch_ncbi_genomeinfo = ch_remote_genome_sources
                .splitCsv(skip: 1, header: true, sep: '\t')
                .map { row ->
                    [
                        accno: row["#assembly_accession"],
                        genome_fna: "${row.ftp_path}/${row.ftp_path - ~/\/$/ - ~/.*\//}_genomic.fna.gz",
                        genome_gff: ""
                    ]
                }

        ch_sample_sigs = channel.empty()
        if ( index_list || ! skip_sourmash ) {
            SAMPLE_SKETCH(ch_sample_reads)
            ch_sample_sigs = SAMPLE_SKETCH.out.signatures
        }

        // Skip sketching and indexing of user-provided genomes if skip_sourmash is set
        ch_joint_user_genomes = ch_user_genomeinfo   // Will be set to selected genomes if sourmash is _not_ skipped, since sourmash will then be used to select matching genomes
        if ( ! skip_sourmash ) {
            GENOME_SKETCH(ch_user_genomeinfo.map { it -> [ [ id: it.accno ], it.genome_fna ] })

            ch_genome_sigs = GENOME_SKETCH.out.signatures
                .collect { _meta, sig -> sig }
                .map { sigs -> [ [ id: 'local-genomes' ], sigs ] }

            GENOME_INDEX(ch_genome_sigs, ksize)

            GATHER_USER_GENOMES(ch_sample_sigs, ch_genome_sigs, true, true, true, true)

            // Collect matching user genomes
            ch_joint_user_genomes = ch_user_genomeinfo
                .map { genome -> [ [ genome.accno ], genome ] }
                .join(
                    GATHER_USER_GENOMES.out.result
                        .map { _meta, csv -> csv }
                        .splitCsv( sep: ',', header: true, quote: '"')
                        .map { genome -> [ [ genome.name ], [ accno: genome.name ] ] }
                        .unique()
                )
                .map { g -> g[1] }
        }

        // Make a cartesian product of samples and user genomes to serve as sample-specific user genomes, i.e. user genomes are always provided to all samples
        ch_sample_user_genomes = ch_sample_reads
            .map { sample -> [ id: sample[0].id ] }
            .combine(ch_joint_user_genomes)

        // Populate both the joint and sample-filtered return channels with all matching user genomes
        ch_joint_filtered_genomes  = ch_joint_user_genomes
        ch_sample_filtered_genomes = ch_joint_user_genomes

        // Call Sourmash with indices for remote genomes if present
        if ( index_list ) {
            // To make sure that all combinations of sample signatures and indexes are gathered below,
            // combine the two channels.
            // (In theory, this should not be required as the command supposedly can take multiple samples
            // and multiple indexes. In our experience, this does not return the full set of hits however.)
            def i = -1
            ch_gather = ch_sample_sigs
                .combine(
                    ch_indexes.map { index -> i += 1; [ [ id: String.format("remoteidx_%02d", i) ], index ] }
                )

            // Do only the remote-genome index gathering here
            GATHER_REMOTE_GENOMES(
                ch_gather.map { it -> [ it[0], it[1] ] },
                ch_gather.map { it -> [ it[2], it[3] ] },
                true, true, true, true
            )

            // The genomes that were selected by sourmash can either be local genomes provided by the user
            // with --genomeinfo, or genomes we need to fetch from NCBI

            // 1. Find the remote genomes that were selected
            // 1.2 Sample-specific set -- "sample"
            ch_sample_remote_genomes = GATHER_REMOTE_GENOMES.out.result
                .splitCsv( sep: ',', header: true, quote: '"')
                // Strip everything except accession number from NCBI-like names
                .map { meta, genome ->
                    def matcher = ( genome.name =~ /^(GC[A-Z]_[0-9]+\.[0-9]+)/ )
                    return matcher ?
                        [ id: meta.id, accno: matcher[0][0] ] :
                        [ id: meta.id, accno: genome.name ]
                }

            // 1.2 Total set -- "joint"
            ch_joint_remote_genomes = ch_sample_remote_genomes
                .map { g -> [ accno: g.accno ] }
                .unique()

            // 2. Fetch NCBI genomes
            WGET_GENOME(
                ch_joint_remote_genomes
                    .map { genome -> [ [ genome.accno ] ] }
                    .join(
                        ch_ncbi_genomeinfo
                            .map { genome -> [ [ genome.accno ], genome ] }
                    )
                    .map { genome -> [ [ id: genome[1].accno ], genome[1].genome_fna ] }
            )

            // 3. Mix the local and the newly fetched NCBI genomes
            // 3.1 Sample specific sets
            ch_sample_filtered_genomes = ch_sample_remote_genomes
                .map { g -> [ [ id: g.accno ], g ] }
                .join(WGET_GENOME.out.outfile)
                .map { g -> [ [  id: g[1].id ], [ accno: g[1].accno, g_fna: g[2] ] ] }
                .mix(ch_sample_user_genomes)

            // 3.2 Total set -- "joint"
            ch_joint_filtered_genomes = ch_joint_user_genomes
                .mix(
                    WGET_GENOME.out.outfile
                        .map { genome -> [ accno: genome[0].id, genome_fna: genome[1] ] }
                )
        }

    emit:
        joint_filtered_genomes  = ch_joint_filtered_genomes
        sample_filtered_genomes = ch_sample_filtered_genomes
}
