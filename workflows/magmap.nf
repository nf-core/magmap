/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAKTA                                   } from '../subworkflows/local/bakta'
include { BAM_SORT_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { BBMAP_ALIGN                            } from '../modules/nf-core/bbmap/align'
include { BBMAP_BBDUK                            } from '../modules/nf-core/bbmap/bbduk'
include { CAT_FASTQ            	                 } from '../modules/nf-core/cat/fastq'
include { CAT_MANY as CAT_GFF                    } from '../modules/local/cat/many'
include { GENOMES2ORFS                           } from '../modules/local/genomes2orfs'
include { CATPROKKATSVS        	                 } from '../modules/local/catprokkatsvs'
include { CHECK_DUPLICATES                       } from '../modules/local/check_duplicates'
include { COLLECT_FEATURECOUNTS                  } from '../modules/local/collect/featurecounts'
include { COLLECT_GENOMESELECTION                } from '../modules/local/collect/genomeselection'
include { TIDYVERSE_JOINFEATURECOUNTSACCNO       } from '../modules/local/tidyverse/joinfeaturecountsaccno'
include { CREATE_BBMAP_INDEX                     } from '../subworkflows/local/create_bbmap_index'
include { CUSTOM_COLLECTSTATS                    } from '../modules/nf-core/custom/collectstats/main'
include { DUCKDB_TABLE2PARQUET                   } from '../modules/nf-core/duckdb/table2parquet'
include { FASTQC                                 } from '../modules/nf-core/fastqc'
include { FASTQC_TRIMGALORE                      } from '../subworkflows/local/fastqc_trimgalore'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { MULTIQC                                } from '../modules/nf-core/multiqc'
include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { PIPELINE_COMPLETION                    } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PIPELINE_INITIALISATION                } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PROKKA                                 } from '../modules/nf-core/prokka'
include { PROKKA_VERSION                         } from '../modules/local/prokka_version'
include { PROKKAGFF2TSV                          } from '../modules/local/prokkagff2tsv'
include { RENAME_CONTIGS                         } from '../modules/local/rename_contigs'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SOURMASH                               } from '../subworkflows/local/sourmash'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts'
include { TIDYVERSE_JOINMETADATA                 } from '../modules/local/tidyverse/joinmetadata/'
include { TIDYVERSE_SELECTANNOTATOR              } from '../modules/local/tidyverse/selectannotator/'
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
    species_preference          //  string: 'all' to select all genomes for a species or 'local', 'completeness' or 'gtdb' to prefer one according to different criteria
    annotator                   //  string: 'prokka', 'bakta_supported_only' or 'bakta_all' -- which tool(s) to annotate genomes lacking a GFF with
    skip_sourmash               // boolean: run Sourmash or not
    sourmash_ksize              // integer
    ch_features                 // channel: list of feature types to call
    skip_fastqc                 // boolean
    skip_qc                     // boolean
    skip_trimming               // boolean
    save_parquet                // boolean: also write summary tables as Parquet
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir                      //  string: output directory

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

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
        species_preference,
        ch_gtdb_metadata,
        ch_gtdbtk_metadata,
        ch_checkm_metadata,
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
    // MODULE: Prokka/Bakta - get gff for all genomes that lack it
    //

    // Find genomes without gff file, and split them between Prokka and Bakta according to
    // --annotator and, for 'bakta_supported_only', each genome's GTDB domain classification
    ch_no_gff = ch_genomes
        .filter { g -> ! g.genome_gff }
        .map { g -> [ [ id: g.accno ], g.genome_fna ] }

    if ( annotator == 'prokka' ) {
        ch_no_gff_prokka = ch_no_gff
        ch_no_gff_bakta  = channel.empty()
    } else {
        TIDYVERSE_SELECTANNOTATOR(
            ch_no_gff.collectFile(name: 'no_gff_accessions.txt', newLine: true) { meta, _fna -> meta.id },
            TIDYVERSE_JOINMETADATA.out.genome_metadata,
            annotator
        )

        // Warn (once) about genomes bakta_supported_only couldn't classify by domain and
        // therefore routed to Prokka rather than Bakta
        TIDYVERSE_SELECTANNOTATOR.out.unclassified_accessions
            .splitText() { it.trim() }
            .filter { it }
            .collect()
            .subscribe { accnos ->
                if ( accnos ) {
                    log.warn "--annotator ${annotator}: could not determine a GTDB domain for ${accnos.size()} genome(s) lacking a GFF (${accnos.join(', ')}) -- routed to Prokka instead of Bakta. Provide --gtdb_metadata/--gtdbtk_metadata covering these genomes to use Bakta for them."
                }
            }

        ch_no_gff_keyed = ch_no_gff.map { meta, fna -> [ meta.id, meta, fna ] }

        ch_no_gff_bakta = ch_no_gff_keyed
            .join(
                TIDYVERSE_SELECTANNOTATOR.out.bakta_accessions
                    .splitText() { it.trim() }
                    .filter { it }
                    .map { accno -> [ accno, true ] }
            )
            .map { _accno, meta, fna, _kept -> [ meta, fna ] }

        ch_no_gff_prokka = ch_no_gff_keyed
            .join(
                TIDYVERSE_SELECTANNOTATOR.out.prokka_accessions
                    .splitText() { it.trim() }
                    .filter { it }
                    .map { accno -> [ accno, true ] }
            )
            .map { _accno, meta, fna, _kept -> [ meta, fna ] }
    }

    PROKKA(ch_no_gff_prokka, [], [])
    ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.log.collect{ _meta, log -> log })
    if ( annotator != 'bakta_all' ) {
        PROKKA_VERSION()
    }

    ch_bakta_fna = channel.empty()
    ch_bakta_gff = channel.empty()
    if ( annotator != 'prokka' ) {
        BAKTA(ch_no_gff_bakta)
        ch_bakta_fna = BAKTA.out.fna
        ch_bakta_gff = BAKTA.out.gff
        ch_multiqc_files = ch_multiqc_files.mix(BAKTA.out.txt.collect{ _meta, txt -> txt })
    }

    PROKKAGFF2TSV(
        ch_genomes.filter { g -> g.genome_gff }.map { g -> [ [ id: g.accno ], g.genome_gff ] }
            .mix(ch_bakta_gff)
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

    // Mix genome entries that were not sent to Prokka/Bakta with those that were
    ch_collected_genomes = ch_genomes
        .filter { g -> g.genome_gff }
        .mix(
            PROKKA.out.fna
                .join(PROKKA.out.gff)
                .map { meta, fna, gff -> [ accno: meta.id  , genome_fna: fna, genome_gff: gff ] }
        )
        .mix(
            ch_bakta_fna
                .join(ch_bakta_gff)
                .map { meta, fna, gff -> [ accno: meta.id  , genome_fna: fna, genome_gff: gff ] }
        )

    if ( genomeset_mode == 'joint' ) {
        //
        // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
        //
        CREATE_BBMAP_INDEX(
            ch_collected_genomes
                .map { it -> [ it.accno, it.genome_fna ] }
                .toList()
                .map { pairs -> [ [ id: 'all' ], pairs.collect { it[0] }, pairs.collect { it[1] } ] }
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
            .map { g -> [ [ id: g[1].id ], g[2].accno, g[2].genome_fna ] }
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
    // Publish the genome accessions that went into each BBMap index -- one file for the
    // whole run in 'joint' mode, one per sample in 'sample' mode. Reused as the input to
    // COLLECT_GENOMESELECTION below.
    //
    ch_genome_accnos_files = CREATE_BBMAP_INDEX.out.genome_accnos
        .collectFile(storeDir: "${outdir}/bbmap") { meta, accnos ->
            [ "${meta.id}.genomes.txt", accnos.sort().join('\n') + '\n' ]
        }

    //
    // MODULE: Summarize local vs remote genome selection for the MultiQC report -- one
    // row per sample when genomeset_mode is 'sample', a single row for the whole run
    // otherwise.
    //
    COLLECT_GENOMESELECTION(
        ch_genome_accnos_files.collect().map { files -> [ [ id: 'magmap' ], files ] },
        SOURMASH.out.local_accessions
    )
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_GENOMESELECTION.out.multiqc.collect { _meta, tsv -> tsv })

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
    ch_fasta_fai = CREATE_BBMAP_INDEX.out.genome_fnas
        .map { meta, fna -> [ meta, fna, [] ] }

    BAM_SORT_STATS_SAMTOOLS(BBMAP_ALIGN.out.bam, ch_fasta_fai)

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
    ch_multiqc_files = ch_multiqc_files.mix(
        FEATURECOUNTS.out.summary
            .map { meta, summary ->
                def content = summary.text.replaceAll(/\S+\.sorted\.bam/, "${meta.id}.${meta.feature}")
                [ "${meta.id}.${meta.feature}.featureCounts.tsv.summary", content ]
            }
            .collectFile { name, content -> [ name, content ] }
            .collect()
        )

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

    COLLECT_FEATURECOUNTS(ch_collect_featurecounts)

    // COLLECT_FEATURECOUNTS itself is kept generic (no genome-accession lookup), so this
    // pipeline-specific join is a separate step -- see nf-core/magmap#237, where this
    // module is being split out into a shared nf-core/modules component and this join
    // stays local, since attaching accno is specific to a genome-collection pipeline.
    // .first() converts the single GENOMES2ORFS emission to a value channel so it's
    // reused for every one of COLLECT_FEATURECOUNTS's per-feature-type emissions --
    // without it, a queue channel with only one item would only pair with the first
    // emission before closing, silently starving the rest.
    TIDYVERSE_JOINFEATURECOUNTSACCNO(
        COLLECT_FEATURECOUNTS.out.counts,
        GENOMES2ORFS.out.genomes2orfs.map { _m, g2orfs -> g2orfs }.first()
    )
    ch_fcs_for_stats      = TIDYVERSE_JOINFEATURECOUNTSACCNO.out.counts.collect { _meta, tsv -> tsv }.map { it -> [ it ] }
    ch_collect_stats      = ch_collect_stats.combine(ch_fcs_for_stats)

    //
    // Collect statistics from the pipeline
    //
    CUSTOM_COLLECTSTATS(ch_collect_stats.map { s -> s + [[]] }) // The last [[]] is to create a value for the `mergetab` that we have in metatdenovo (which shares the swf)

    //
    // MODULE: Also write the summary tables as Parquet
    //
    if ( save_parquet ) {
        DUCKDB_TABLE2PARQUET(
            TIDYVERSE_JOINMETADATA.out.genome_metadata
                .mix(COLLECT_GENOMESELECTION.out.full_table.map { _meta, tsv -> tsv })
                .mix(GENOMES2ORFS.out.genomes2orfs.map { _meta, tsv -> tsv })
                .mix(CATPROKKATSVS.out.tsv.map { _meta, tsv -> tsv })
                .mix(TIDYVERSE_JOINFEATURECOUNTSACCNO.out.counts.map { _meta, tsv -> tsv })
                .mix(CUSTOM_COLLECTSTATS.out.overall_stats.map { _meta, tsv -> tsv })
                .map { tsv -> [ [ id: tsv.name.replaceAll(/\.tsv(\.gz)?$/, '') ], tsv ] }
        )
    }

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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
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
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'magmap'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
