/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_SORT_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { BBMAP_ALIGN                            } from '../modules/nf-core/bbmap/align'
include { BBMAP_BBDUK                            } from '../modules/nf-core/bbmap/bbduk'
include { CAT_FASTQ            	                 } from '../modules/nf-core/cat/fastq'
include { CATPROKKATSVS        	                 } from '../modules/local/catprokkatsvs'
include { CONCATENATE_GFFS                       } from '../subworkflows/local/concatenate_gffs'
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
    ch_versions = ch_versions.mix(CHECK_DUPLICATES.out.versions)

    ch_duplicates = CHECK_DUPLICATES.out.duplicate_genomes
        .flatMap { it.tokenize('\n') }
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
    ch_versions = ch_versions.mix(RENAME_CONTIGS.out.versions)

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
        .branch { meta, reads ->
            cat: reads.size() >= 2
            skip_cat: true
        }

    //
    // MODULE: Run FastQC on the raw reads
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_short_reads_forcat.cat.map { meta, reads -> [meta, reads.flatten()] })
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    // Ensure we don't have nests of nests so that structure is in form expected for assembly
    ch_short_reads_catskipped = ch_short_reads_forcat.skip_cat.map { meta, reads ->
        def new_reads = meta.single_end ? reads[0] : reads.flatten()
            [meta, new_reads]
    }

    // Combine single run and multi-run-merged data
    ch_short_reads = channel.empty()
    ch_short_reads = CAT_FASTQ.out.reads.mix(ch_short_reads_catskipped)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_short_reads,
        skip_fastqc || skip_qc,
        skip_trimming
    )

    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    ch_collect_stats = ch_short_reads
        .collect { meta, fasta -> meta }
        .map { reads -> [ [ id: 'magmap' ], reads ] }

    if ( skip_trimming ) {
        ch_collect_stats = ch_collect_stats
            .map { meta, samples -> [ meta, samples, [] ] }

    } else {
        ch_collect_stats = ch_collect_stats
            .combine(
                FASTQC_TRIMGALORE.out.trim_log
                    .collect { meta, report ->
                        if ( report in List ) {
                            report[0]
                        } else {
                            report
                        }
                    }
                    .map { [ it ] }
            )
    }

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via --sequence_filter
    //
    if ( sequence_filter ) {
        BBMAP_BBDUK(FASTQC_TRIMGALORE.out.reads, sequence_filter)
        ch_versions   = ch_versions.mix(BBMAP_BBDUK.out.versions)

        ch_clean_reads = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { it[1] }.map { [ it ] }
        ch_collect_stats = ch_collect_stats
            .combine(ch_bbduk_logs)
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{ meta, log -> log })
    } else {
        ch_clean_reads = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = channel.empty()
        ch_collect_stats = ch_collect_stats
            .map { [ it[0], it[1], it[2], [] ] }
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
    ch_versions = ch_versions.mix(SOURMASH.out.versions)
    ch_genomes = SOURMASH.out.filtered_genomes

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
    ch_versions = ch_versions.mix(TIDYVERSE_JOINMETADATA.out.versions.first())

    //
    // MODULE: Prokka - get gff for all genomes that lack it
    //

    // Find genomes without gff file, and pass them to Prokka
    ch_no_gff = ch_genomes
        .filter { g -> ! g.genome_gff }
        .map { g -> [ [ id: g.accno ], g.genome_fna ] }

    PROKKA(ch_no_gff, [], [])
    ch_versions = ch_versions.mix(PROKKA.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.log.collect{ meta, log -> log })

    PROKKAGFF2TSV(
        ch_genomes.filter { g -> g.genome_gff }.map { g -> [ [ id: g.accno ], g.genome_gff ] }
    )
    ch_versions = ch_versions.mix(PROKKAGFF2TSV.out.versions)

    CATPROKKATSVS(
        PROKKA.out.tsv
            .map { t -> t[1] }
            .mix(
                PROKKAGFF2TSV.out.tsv.map { t -> t[1] }
            )
            .collect()
            .map { t -> [ [ id: 'magmap' ], t ] }
    )
    ch_versions = ch_versions.mix(CATPROKKATSVS.out.versions)

    // Mix genome entries that were not sent to Prokka with those that were not
    ch_collected_genomes = ch_genomes
        .filter { g -> g.genome_gff }
        .mix(
            PROKKA.out.fna
                .join(PROKKA.out.gff)
                .map { meta, fna, gff -> [ accno: meta.id  , genome_fna: fna, genome_gff: gff ] }
        )

    //
    // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
    //
    CREATE_BBMAP_INDEX(ch_collected_genomes.map { it.genome_fna })
    ch_versions = ch_versions.mix(CREATE_BBMAP_INDEX.out.versions)

    //
    // SUBWORKFLOW: Concatenate gff files
    //
    CONCATENATE_GFFS(ch_collected_genomes.map { it.genome_gff })
    ch_versions = ch_versions.mix(CONCATENATE_GFFS.out.versions)

    //
    // BBMAP ALIGN. Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, CREATE_BBMAP_INDEX.out.index )
    ch_multiqc_files = ch_multiqc_files.mix(BBMAP_ALIGN.out.log.collect{ meta, log -> log })
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    //
    // SUBWORKFLOW: sort bam file and produce statistics
    //
    BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, CREATE_BBMAP_INDEX.out.genome_fnas)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    ch_stage_counts = BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(CONCATENATE_GFFS.out.gff.map { it[1] })

    ch_collect_stats = ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it[1]}.map { [ it ] })

    //
    // MODULE: FeatureCounts
    //
    ch_featurecounts = ch_stage_counts
        .combine(ch_features)
        .map { meta, bam, gff, feature ->
            [ meta + [feature: feature], bam, gff ]
        }

    FEATURECOUNTS(ch_featurecounts)
    ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.collect{ meta, log -> log })

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    ch_collect_featurecounts = FEATURECOUNTS.out.counts
        .map { meta, file -> [ meta.feature, [meta, file] ] }
        .groupTuple()
        .map { feature, data ->
            def metas = data.collect { it[0] }
            def files = data.collect { it[1] }
            [ metas[0] + [feature: feature], files ]
        }
        .map { meta, data ->
            [ [id: meta.feature ], data ]
        }

    COLLECT_FEATURECOUNTS(ch_collect_featurecounts)
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)

    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { meta, tsv -> tsv }.map { [ it ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { meta, tsv -> tsv }
    ch_collect_stats      = ch_collect_stats.combine(ch_fcs_for_stats)

    //
    // Collect statistics from the pipeline
    //
    COLLECT_STATS(ch_collect_stats.map { s -> s + [[]] }) // The last [[]] is to create a value for the `mergetab` that we have in metatdenovo (which shares the swf)
    ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

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

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
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

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
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
