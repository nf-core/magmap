/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_SORT_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { BBMAP_ALIGN                            } from '../modules/nf-core/bbmap/align'
include { BBMAP_BBDUK                            } from '../modules/nf-core/bbmap/bbduk'
include { CAT_FASTQ            	                 } from '../modules/nf-core/cat/fastq'
include { CAT_MANY as CAT_GFF                    } from '../modules/local/cat/many'
include { GENOMES2ORFS                           } from '../modules/local/genomes2orfs'
include { CATPROKKATSVS        	                 } from '../modules/local/catprokkatsvs'
include { CHECK_DUPLICATES                       } from '../modules/local/check_duplicates'
include { COLLECT_FEATURECOUNTS                  } from '../modules/local/collect/featurecounts'
include { COLLECT_STATS                          } from '../modules/local/collect/stats'
include { CREATE_BBMAP_INDEX                     } from '../subworkflows/local/create_bbmap_index'
include { FASTQC                                 } from '../modules/nf-core/fastqc'
include { FASTQC_TRIMGALORE                      } from '../subworkflows/local/fastqc_trimgalore'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { MULTIQC                                } from '../modules/nf-core/multiqc'
include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { PIPELINE_COMPLETION                    } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PIPELINE_INITIALISATION                } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PROKKA                                 } from '../modules/nf-core/prokka'
include { PROKKAGFF2TSV                          } from '../modules/local/prokkagff2tsv'
include { RENAME_CONTIGS                         } from '../modules/local/rename_contigs'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SOURMASH                               } from '../subworkflows/local/sourmash'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts'
include { TIDYVERSE_JOINMETADATA                 } from '../modules/local/tidyverse/joinmetadata/'
include { validateInputSamplesheet               } from '../subworkflows/local/utils_nfcore_magmap_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAGMAP {

    take:
    ch_samplesheet              // channel: samplesheet read in from --input
    ch_genomeinfo               // channel: genome information sheet read in from --genomeinfo
    ch_remote_genome_sources    // channel: paths to NCBI-style genome summary files
    ch_indexes                  // channel: user-provided Sourmash indexes
    index_list                  //  string: value of the indexes params, used for if clauses in the SOURMASH subworkflow
    sequence_filter             //  string: fasta file for BBDuk
    ch_gtdb_metadata            // channel: GTDB metadata files
    ch_gtdbtk_metadata          // channel: GTDB-Tk metadata files
    ch_checkm_metadata          // channel: CheckM/CheckM2 metadata files
    genomeset_mode              //  string: Either 'joint' for mapping samples against all genomes, or 'sample' to map to sample-specific sets
    skip_sourmash               // boolean: run Sourmash or not
    sourmash_ksize              // integer
    ch_features                 // channel: list of feature types to call
    skip_fastqc                 // boolean
    skip_qc                     // boolean
    skip_trimming               // boolean
    outdir                      //  string: output directory

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Check presence of duplicates contigs in the local genome collection
    //
    ch_check_duplicates = ch_genomeinfo
        .collect { g -> g.genome_fna }
        .map { g -> [ [ id: "local-genomes" ], g ] }

    CHECK_DUPLICATES(ch_check_duplicates)

    ch_duplicates = CHECK_DUPLICATES.out.duplicate_genomes
        .flatMap { it -> it.tokenize('\n') }
        .map { fname -> [ fname.replaceAll(/.*\//, ''), true ] }
    ch_genomes_pre_renaming = ch_genomeinfo
        .map { row -> [ row.genome_fna.getName(), row ] }
        .join(ch_duplicates, remainder: true)
        .map { row -> [ row[1], row[2] ] }
        .branch { row ->
            needs_renaming: row[1]
                return row[0]
            names_ok:       true
                return row[0]
        }

    RENAME_CONTIGS(
        ch_genomes_pre_renaming.needs_renaming
            .map { g -> [ [ id: g.accno ], g.genome_fna ] }
    )

    ch_genomes_post_renaming = RENAME_CONTIGS.out.renamed_contigs
        .map { g -> [ accno: g[0].id, genome_fna: g[1], genome_gff: [] ] }
        .mix(ch_genomes_pre_renaming.names_ok)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_short_reads_forcat = ch_samplesheet
        .map { meta, reads ->
            def meta_new = meta - meta.subMap('run')
            [meta_new, reads]
        }
        .groupTuple()
        .branch { _meta, reads ->
            cat: reads.size() >= 2
            skip_cat: true
        }

    //
    // MODULE: Run FastQC on the raw reads
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it -> it[1] })

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_short_reads_forcat.cat.map { meta, reads -> [meta, reads.flatten()] })

    // Ensure we don't have nests of nests so that structure is in form expected for assembly
    ch_short_reads_catskipped = ch_short_reads_forcat.skip_cat.map { meta, reads ->
        def new_reads = meta.single_end ? reads[0] : reads.flatten()
            [meta, new_reads]
    }

    // Combine single run and multi-run-merged data
    ch_short_reads = channel.empty()
    ch_short_reads = CAT_FASTQ.out.reads.mix(ch_short_reads_catskipped)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_short_reads,
        skip_fastqc || skip_qc,
        skip_trimming
    )

    ch_collect_stats = ch_short_reads
        .collect { meta, _fasta -> meta }
        .map { reads -> [ [ id: 'magmap' ], reads ] }

    if ( skip_trimming ) {
        ch_collect_stats = ch_collect_stats
            .map { meta, samples -> [ meta, samples, [] ] }

    } else {
        ch_collect_stats = ch_collect_stats
            .combine(
                FASTQC_TRIMGALORE.out.trim_log
                    .collect { _meta, report ->
                        if ( report in List ) {
                            report[0]
                        } else {
                            report
                        }
                    }
                    .map { it -> [ it ] }
            )
    }

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via --sequence_filter
    //
    if ( sequence_filter ) {
        BBMAP_BBDUK(FASTQC_TRIMGALORE.out.reads, sequence_filter)

        ch_clean_reads = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { it -> it[1] }.map { it -> [ it ] }
        ch_collect_stats = ch_collect_stats
            .combine(ch_bbduk_logs)
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{ _meta, log -> log })
    } else {
        ch_clean_reads = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = channel.empty()
        ch_collect_stats = ch_collect_stats
            .map { it -> [ it[0], it[1], it[2], [] ] }
    }

    //
    // SUBWORKFLOW: Use SOURMASH on sample reads and genomes to reduce the number of the latter
    //
    SOURMASH(
        ch_clean_reads,
        ch_indexes,
        index_list,
        ch_genomes_post_renaming,
        ch_remote_genome_sources,
        sourmash_ksize,
        skip_sourmash
    )
    ch_genomes = SOURMASH.out.joint_filtered_genomes

    //
    // MODULE: Join and filter genome metadata
    //
    TIDYVERSE_JOINMETADATA(
        ch_genomes
            .collectFile(
                name: 'selected_genomes.tsv',
                newLine: true
            ) { genome_record -> genome_record.accno },
        ch_gtdb_metadata.collect().ifEmpty([]),
        ch_gtdbtk_metadata.collect().ifEmpty([]),
        ch_checkm_metadata.collect().ifEmpty([])
    )

    //
    // MODULE: Prokka - get gff for all genomes that lack it
    //

    // Find genomes without gff file, and pass them to Prokka
    ch_no_gff = ch_genomes
        .filter { g -> ! g.genome_gff }
        .map { g -> [ [ id: g.accno ], g.genome_fna ] }

    PROKKA(ch_no_gff, [], [])
    ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.log.collect{ _meta, log -> log })

    PROKKAGFF2TSV(
        ch_genomes.filter { g -> g.genome_gff }.map { g -> [ [ id: g.accno ], g.genome_gff ] }
    )

    CATPROKKATSVS(
        PROKKA.out.tsv
            .map { t -> t[1] }
            .mix(
                PROKKAGFF2TSV.out.tsv.map { t -> t[1] }
            )
            .collect()
            .map { t -> [ [ id: 'magmap' ], t ] }
    )

    // Mix genome entries that were not sent to Prokka with those that were not
    ch_collected_genomes = ch_genomes
        .filter { g -> g.genome_gff }
        .mix(
            PROKKA.out.fna
                .join(PROKKA.out.gff)
                .map { meta, fna, gff -> [ accno: meta.id  , genome_fna: fna, genome_gff: gff ] }
        )

    if ( genomeset_mode == 'joint' ) {
        //
        // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
        //
        CREATE_BBMAP_INDEX(
            ch_collected_genomes
                .collect { it -> it.genome_fna }
                .map { it -> [ [ id: 'all' ], it ] }
        )

        //
        // BBMAP ALIGN. Call BBMap with the index once per sample
        //
        BBMAP_ALIGN(ch_clean_reads, CREATE_BBMAP_INDEX.out.index.map { index -> index[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_ALIGN.out.log.collect{ _meta, log -> log })
    } else if ( genomeset_mode == 'sample' ) {
        ch_fnas_to_index = SOURMASH.out.sample_filtered_genomes
            .map { g -> [ [ accno: g[1].accno ], [ id: g[0].id, accno: g[1].accno ] ] }
            .combine(ch_collected_genomes.map { g -> [ [ accno: g.accno ], g ] }, by: 0)
            .map { g -> [ [ id: g[1].id ], g[2].genome_fna ] }
            .groupTuple()

        //
        // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
        //
        CREATE_BBMAP_INDEX(ch_fnas_to_index)

        //
        // BBMAP ALIGN. Call BBMap with the index once per sample
        //

        // Make sure the correct index is sent with each sample
        ch_reads_and_indices = ch_clean_reads
            .map { r -> [ [ id: r[0].id ], r[0], r[1] ] }
            .join(CREATE_BBMAP_INDEX.out.index)

        BBMAP_ALIGN(
            ch_reads_and_indices.map { ri -> [ ri[1], ri[2] ] },
            ch_reads_and_indices.map { ri -> ri[3] }
        )
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_ALIGN.out.log.collect { _meta, log -> log })
    }

    //
    // MODULE: Concatenate gff files
    //
    CAT_GFF([id: 'genomes'], ch_collected_genomes.map { genome -> genome.genome_gff }.collect())

    //
    // MODULE: Create an index file from genome accnos to feature prefixes
    //
    GENOMES2ORFS(ch_collected_genomes.map { genome -> genome.genome_gff }.collect().map { genomes -> [ [ id: 'genomes' ], genomes ] })

    //
    // SUBWORKFLOW: sort bam file and produce statistics
    //
    BAM_SORT_STATS_SAMTOOLS(BBMAP_ALIGN.out.bam, CREATE_BBMAP_INDEX.out.genome_fnas)

    ch_stage_counts = BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(CAT_GFF.out.concatenated_files.map { it -> it[1] })

    ch_collect_stats = ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it -> it[1]}.map { it -> [ it ] })

    //
    // MODULE: FeatureCounts
    //
    ch_featurecounts = ch_stage_counts
        .combine(ch_features)
        .map { meta, bam, gff, feature ->
            [ meta + [feature: feature], bam, gff ]
        }

    FEATURECOUNTS(ch_featurecounts)
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.collect{ _meta, log -> log })

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    ch_collect_featurecounts = FEATURECOUNTS.out.counts
        .map { meta, file -> [ meta.feature, [meta, file] ] }
        .groupTuple()
        .map { feature, data ->
            def metas = data.collect { it -> it[0] }
            def files = data.collect { it -> it[1] }
            [ metas[0] + [feature: feature], files ]
        }
        .map { meta, data ->
            [ [id: meta.feature ], data ]
        }

    COLLECT_FEATURECOUNTS(ch_collect_featurecounts, GENOMES2ORFS.out.genomes2orfs.map { _m, g2orfs -> g2orfs })

    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { _meta, tsv -> tsv }.map { it -> [ it ] }
    //ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { _meta, tsv -> tsv }
    ch_collect_stats      = ch_collect_stats.combine(ch_fcs_for_stats)

    //
    // Collect statistics from the pipeline
    //
    COLLECT_STATS(ch_collect_stats.map { s -> s + [[]] }) // The last [[]] is to create a value for the `mergetab` that we have in metatdenovo (which shares the swf)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'magmap_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    //ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
