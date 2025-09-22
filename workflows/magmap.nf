/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_SORT_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BBMAP_ALIGN                            } from '../modules/nf-core/bbmap/align/main'
include { BBMAP_BBDUK                            } from '../modules/nf-core/bbmap/bbduk/main'
include { CAT_GFFS                               } from '../subworkflows/local/concatenate_gff'
include { CAT_FASTQ            	                 } from '../modules/nf-core/cat/fastq/main'
include { CHECK_DUPLICATES                       } from '../modules/local/check_duplicates'
include { COLLECT_FEATURECOUNTS                  } from '../modules/local/collect_featurecounts'
include { COLLECT_STATS                          } from '../modules/local/collect_stats'
include { CREATE_BBMAP_INDEX                     } from '../subworkflows/local/create_bbmap_index'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { FASTQC_TRIMGALORE                      } from '../subworkflows/local/fastqc_trimgalore'
include { FILTER_GENOMES                         } from '../modules/local/filter_genomes'
include { FORMAT_KRONA                           } from '../modules/local/format_krona/main'
include { KRAKEN2_KRAKEN2                        } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA              } from '../modules/nf-core/krakentools/kreport2krona/main'
include { KRAKEN2_DOWNLOAD_DB                    } from '../modules/local/kraken2/download/main'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { PIGZ_COMPRESS as PIGZ_PE_READS_FWD     } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_PE_READS_REV     } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_SE_READS         } from '../modules/nf-core/pigz/compress/main'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { PIGZ_UNCOMPRESS as GUNZIP_CONTIGS      } from '../modules/nf-core/pigz/uncompress/main'
include { PIPELINE_COMPLETION                    } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PIPELINE_INITIALISATION                } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PROKKA                                 } from '../modules/nf-core/prokka/main'
include { RENAME_CONTIGS                         } from '../modules/local/rename_contigs'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SOURMASH                               } from '../subworkflows/local/sourmash'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'
include { TAXBURST                               } from '../modules/local/taxburst/main'
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
    sequence_filter             //  string: fasta file for BBDuk
    ch_gtdb_metadata            // channel: GTDB metadata files
    ch_gtdbtk_metadata          // channel: GTDB-Tk metadata files
    ch_checkm_metadata          // channel: CheckM/CheckM2 metadata files
    skip_kraken2                // boolean: run Kraken2 or not
    kraken2_db                  // string: path to Kraken2 database
    kraken2_db_url              // string: URL to download Kraken2 database
    sourmash                    // boolean: run Sourmash or not
    sourmash_ksize              // integer
    sourmash_save_unassigned    // boolean
    sourmash_save_matches_sig   // boolean
    sourmash_save_prefetch      // boolean
    sourmash_save_prefetch_csv  // boolean
    ch_features                 // channel: list of feature types to call
    skip_fastqc                 // boolean
    skip_qc                     // boolean
    skip_trimming               // boolean
    outdir                      //  string: output directory

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Check presence of duplicates contigs in the local genome collection
    //
    ch_check_duplicates = ch_genomeinfo
        .map { it.genome_fna }
        .collect()
        .map { [ [id: 'genomes'], it ] }

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

        // Form a fastq channel from the samplesheet channel
    // DL: I'm not sure which parts are still required after nf-schema. The branch { } certainly is needed.
    ch_fastq = ch_samplesheet
        .flatMap { meta, fastq_files ->
            if (fastq_files.size() <= 2) {
                return [[ meta.id, [meta], fastq_files ]]
            } else {
                def pairs = fastq_files.collate(2)
                return [[ meta.id, pairs.collect { meta + [id: "${meta.id}_${pairs.indexOf(it) + 1}"] }, fastq_files ]]
            }
        }
        /** DL: In my testing, this fails as entries come in with single meta, multiple read files when appear for single ends
        .map { id, metas, fastq_files ->
            // Ensure single_end is set correctly in meta
            def updatedMetas = metas
                .collect { it + [single_end: (fastq_files.size() / metas.size() == 1)] }
            return [id, updatedMetas, fastq_files]
        }
        **/
        .map { validateInputSamplesheet(it) }
        .branch {
            meta, fastqs ->
                single  : ( meta.single_end && fastqs.size() == 1 ) || ( ! meta.single_end && fastqs.size == 2 )
                    return [ meta, fastqs ]
                multiple: true
                    return [ meta, fastqs ]
        }

    //
    // MODULE: Concatenate FastQ files from the same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // Gzip unzipped read files
    //
    // We're only doing this for samples having a single row in the sample sheet since those with more than one
    // were gzipped by CAT_FASTQ above.
    //

    // Paired end, forward
    fwd = ch_fastq.single
        .filter { meta, f -> ! meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[0] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_PE_READS_FWD(fwd.unzipped)
    ch_versions      = ch_versions.mix(PIGZ_PE_READS_FWD.out.versions)

    // Paired end, reverse
    rev = ch_fastq.single
        .filter { meta, f -> ! meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[1] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_PE_READS_REV(rev.unzipped)
    ch_versions      = ch_versions.mix(PIGZ_PE_READS_REV.out.versions)

    // Single end
    se = ch_fastq.single
        .filter { meta, f -> meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[0] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_SE_READS(se.unzipped)
    ch_versions      = ch_versions.mix(PIGZ_SE_READS.out.versions)

    // Join the three channels with the originally zipped to form a new ch_fastq of the same structure as the original
    ch_fastq = fwd.zipped.concat(PIGZ_PE_READS_FWD.out.archive)
        .join(rev.zipped.concat(PIGZ_PE_READS_REV.out.archive))
        .map { meta, fwd, rev -> [ meta, [ fwd, rev ] ] }
        .concat(
            se.zipped
                .concat(PIGZ_SE_READS.out.archive)
                .map { meta, fastq -> [ meta, [ fastq ] ] }
        )
        .concat(CAT_FASTQ.out.reads)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)


    ch_collect_stats = ch_fastq
        .collect { 
            meta, fasta -> meta.id 
            }
        .map {
            [ [ id:"magmap" ], it ] 
        }

    if ( params.skip_trimming ) {
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
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, Channel.fromPath(params.sequence_filter).first() )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { meta, log ->  log }.map { [ it ] }
        ch_versions   = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        ch_collect_stats = ch_collect_stats.combine(ch_bbduk_logs)
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{ meta, log -> log })
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = Channel.empty()
        ch_collect_stats = ch_collect_stats
            .map { meta, samples, report -> [ meta, samples, report, [] ] }
    }

    //
    // MODULE: Kraken2
    //
    if ( ! skip_kraken2 ) {
        if (!kraken2_db) {

            db_url = kraken2_db_url

            KRAKEN2_DOWNLOAD_DB(db_url)
            ch_kraken2_db = KRAKEN2_DOWNLOAD_DB.out.db_dir
            ch_versions = ch_versions.mix(KRAKEN2_DOWNLOAD_DB.out.versions)

        } else {
            ch_kraken2_db = Channel.fromPath(params.kraken2_db, checkIfExists: true, type: 'dir')
        }

        KRAKEN2_KRAKEN2(
            ch_clean_reads,
            ch_kraken2_db,
            params.kraken2_save_output,
            params.kraken2_save_reads
        )

        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

        KRAKENTOOLS_KREPORT2KRONA(KRAKEN2_KRAKEN2.out.report)
        ch_versions = ch_versions.mix(KRAKENTOOLS_KREPORT2KRONA.out.versions)

        FORMAT_KRONA(KRAKENTOOLS_KREPORT2KRONA.out.txt)
        ch_versions = ch_versions.mix(FORMAT_KRONA.out.versions)

        TAXBURST(FORMAT_KRONA.out.format_tsv)
        ch_versions = ch_versions.mix(TAXBURST.out.versions)
    }

    //
    // SUBWORKFLOW: Use SOURMASH on sample reads and genomes to reduce the number of the latter
    //
    if ( sourmash ) {
        SOURMASH(
            ch_clean_reads,
            ch_indexes,
            ch_genomes_post_renaming,
            ch_remote_genome_sources,
            sourmash_ksize,
            sourmash_save_unassigned,
            sourmash_save_matches_sig,
            sourmash_save_prefetch,
            sourmash_save_prefetch_csv
        )
        ch_versions = ch_versions.mix(SOURMASH.out.versions)
        ch_genomes = SOURMASH.out.filtered_genomes
    } else {
        ch_genomes = ch_genomes_post_renaming
    }

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

    // First, contigs need to be gunzipped
    ch_no_gff = ch_genomes
        .filter { g -> ! g.genome_gff }
        .map { g -> [ [ id: g.accno ], g.genome_fna ] }
        .branch { g ->
            gzipped: g[1] =~ /\.gz$/
            unzipped: true
        }
    GUNZIP_CONTIGS(ch_no_gff.gzipped)
    ch_versions = ch_versions.mix(GUNZIP_CONTIGS.out.versions)

    PROKKA(GUNZIP_CONTIGS.out.file.mix(ch_no_gff.unzipped), [], [])
    ch_versions = ch_versions.mix(PROKKA.out.versions)

    // PROKKA on the genomes that lack gff
    ch_finished_genomes = ch_genomes
        .filter { g -> g.genome_gff }
        .mix(
            PROKKA.out.gff
                .map{ meta, gff -> [ meta.id  , [ meta.id, gff ] ] }
                .join(ch_no_gff.gzipped.mix(ch_no_gff.unzipped).map { meta, fna -> [ meta.id , [ meta.id, fna ] ] } )
                .map{ meta, gff, fna -> [ accno: gff[0], genome_fna: fna[1], genome_gff: gff[1] ] }
        )

    //
    // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
    //
    CREATE_BBMAP_INDEX ( ch_finished_genomes.map{ it.genome_fna } )
    ch_versions = ch_versions.mix(CREATE_BBMAP_INDEX.out.versions)

    //
    // SUBWORKFLOW: Concatenate gff files
    //
    CAT_GFFS ( ch_finished_genomes.map{ it.genome_gff } )
    ch_versions = ch_versions.mix(CAT_GFFS.out.versions)

    //
    // BBMAP ALIGN. Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, CREATE_BBMAP_INDEX.out.index )
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    //
    // SUBWORKFLOW: sort bam file and produce statistics
    //
    BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, CREATE_BBMAP_INDEX.out.genome_fnas)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    ch_stage_counts = BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(CAT_GFFS.out.gff.map { it[1] })

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

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    ch_collect_featurecounts = FEATURECOUNTS.out.counts
        .map { meta, file -> [meta.feature, [meta, file]] }
        .groupTuple()
        .map { feature, data ->
            def metas = data.collect { it[0] }
            def files = data.collect { it[1] }
            [metas[0] + [feature: feature], files]
        }
        .map { meta, data ->
            [ [id: meta.feature ], data ]
        }

    COLLECT_FEATURECOUNTS ( ch_collect_featurecounts )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { it[1]}.map { [ it ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { it[1]}
    ch_collect_stats = ch_collect_stats
        .combine(ch_fcs_for_stats)

    //
    // Collect statistics from the pipeline
    //
    COLLECT_STATS(ch_collect_stats)
    ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    //
    // Collate and save software versions
    //
    ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'magmap_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
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
