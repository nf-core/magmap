/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLLECT_FEATURECOUNTS                  } from '../modules/local/collect_featurecounts'
include { COLLECT_STATS                          } from '../modules/local/collect_stats'
include { FILTER_GENOMES                         } from '../modules/local/filter_genomes'
include { CHECK_DUPLICATES                       } from '../modules/local/check_duplicates'
include { validateInputSamplesheet               } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { FASTQC_TRIMGALORE                      } from '../subworkflows/local/fastqc_trimgalore'
include { CAT_GFFS                               } from '../subworkflows/local/concatenate_gff'
include { CREATE_BBMAP_INDEX                     } from '../subworkflows/local/create_bbmap_index'
include { SOURMASH                               } from '../subworkflows/local/sourmash'
include { PIPELINE_INITIALISATION                } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PIPELINE_COMPLETION                    } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { BBMAP_BBDUK                            } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_ALIGN                            } from '../modules/nf-core/bbmap/align/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'
include { GUNZIP                                 } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFFS                  } from '../modules/nf-core/gunzip/main'
include { PROKKA                                 } from '../modules/nf-core/prokka/main'
include { CAT_FASTQ            	                 } from '../modules/nf-core/cat/fastq/main'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_SORT_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { UTILS_NEXTFLOW_PIPELINE                } from '../subworkflows/nf-core/utils_nextflow_pipeline/main'
include { UTILS_NFCORE_PIPELINE                  } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAGMAP {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // INPUT: if user provides, populate ch_genomeinfo with a table that provides the genomes to filter with sourmash
    //
    ch_genomeinfo = Channel.empty()
    if ( params.genomeinfo) {
        Channel
            .fromPath( params.genomeinfo )
            .splitCsv( sep: ',', header: true )
            .set { ch_genomeinfo }
    }

    //
    // Check presence of duplicates contigs in the local genome collection
    //
    CHECK_DUPLICATES( ch_genomeinfo.map{ it.genome_fna }.collect() )
    ch_versions = ch_versions.mix(CHECK_DUPLICATES.out.versions)

    //
    // INPUT: genome info from ncbi
    //
    if ( params.ncbi_genome_infos) {
        Channel
            .fromPath( params.ncbi_genome_infos )
            .set { ch_genome_infos }
    }

    //
    // INPUT: if user provides, populate ch_indexes
    //
    ch_indexes = Channel.empty()

    if ( params.indexes) {
        Channel
            .fromPath( params.indexes )
            .set { ch_indexes }
    }

    //
    // INPUT: if user provides, populate ch_metadata
    //
    ch_gtdb_metadata = Channel.empty()

    if ( params.gtdb_metadata) {
        Channel
            .of(params.gtdb_metadata.split(','))
            .map { file(it) }
            .splitCsv( sep: '\t', header: true)
            .map {
                [
                    accno: it.accession - ~/^[A-Z][A-Z]_/,
                    checkm_completeness: it.checkm_completeness,
                    checkm_contamination: it.checkm_contamination,
                    checkm_strain_heterogeneity: it.checkm_strain_heterogeneity,
                    contig_count: it.contig_count,
                    genome_size: it.genome_size,
                    gtdb_genome_representative: it.gtdb_genome_representative,
                    gtdb_representative: it.gtdb_representative,
                    gtdb_taxonomy: it.gtdb_taxonomy
                ]
            }
            .set { ch_gtdb_metadata }
    }

    //
    // INPUT: CheckM metadata
    //
    ch_checkm_metadata = Channel.empty()
    if ( params.checkm_metadata) {
        Channel
            .of(params.checkm_metadata.split(','))
            .map { file(it) }
            .splitCsv( sep: '\t', header: true)
            .map { [ [ it["Bin Id"] ],
                [
                    checkm_completeness: it.Completeness,
                    checkm_contamination: it.Contamination,
                    contig_count: it["# contigs"],
                    checkm_strain_heterogeneity: it["Strain heterogeneity"],
                    genome_size: it["Genome size (bp)"]
                ] ]
            }
            .set { ch_checkm_metadata }
    }

    //
    // INPUT: GTDB-Tk metadata
    //
    ch_gtdbtk_metadata = Channel.empty()
    if ( params.gtdbtk_metadata) {
        Channel
            .of(params.gtdbtk_metadata.split(','))
            .map { file(it) }
            .splitCsv( sep: '\t', header: true)
            .map { [ [it.user_genome],
                [
                    gtdb_genome_representative: it.user_genome,
                    gtdb_representative: "f",
                    gtdb_taxonomy: it.classification
                ] ]
            }
            .set { ch_gtdbtk_metadata }
    }

    //
    // gtdbtk_metadata and checkm_metadata need to be joined
    //
    if ( params.gtdbtk_metadata && params.checkm_metadata && params.gtdb_metadata) {
        ch_gtdbtk_metadata
        .map{ accno, gtdbtk -> [ accno, gtdbtk ]}
        .join( ch_checkm_metadata
            .map{ accno, checkm -> [ accno, checkm ] }
        )
        .map { accno, gtdbtk, checkm ->
                [
                    accno,
                    [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                    ]
                ]
            }
    .set { ch_checkm_gtdb_metadata }

    ch_checkm_gtdb_metadata
            .map {
                accno, gtdbtk_checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 !in it }
            .map{ it[0] }
            .join(
                ch_checkm_gtdb_metadata.map { accno, gtdbtk_checkm -> [ gtdbtk_checkm.accno, gtdbtk_checkm ] }
            )
            .map{ it[1]}
        .set{ ch_gtdbtk_checkm_filtered }

        ch_checkm_gtdb_metadata
            .map {
                accno, gtdbtk_checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 in it }
            .map{ it[0] }
            .join(
                ch_gtdb_metadata.map { gtdb -> [ gtdb.accno, gtdb ] }
            )
            .map{ it[1]}
            .mix(ch_gtdbtk_checkm_filtered)
            .set{ ch_metadata }

    } else if ( !params.gtdbtk_metadata && params.checkm_metadata && params.gtdb_metadata) {
        ch_checkm_metadata
            .map{ accno, checkm -> [
                    accno,
                    [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: "",
                    gtdb_representative: "",
                    gtdb_taxonomy: ""
                    ]
                ]
            }
            .set { ch_checkm_metadata }

        ch_checkm_metadata
            .map {
                accno, checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 !in it }
            .map{ it[0] }
            .join(
                ch_checkm_metadata.map { accno, checkm -> [ checkm.accno, checkm ] }
            )
            .map{ it[1]}
        .set{ ch_checkm_filtered }

        ch_checkm_metadata
            .map {
                accno, checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 in it }
            .map{ it[0] }
            .join(
                ch_gtdb_metadata.map { gtdb -> [ gtdb.accno, gtdb ] }
            )
            .map{ it[1] }
            .mix(ch_checkm_filtered)
        .set{ ch_metadata }

    } else if( params.gtdbtk_metadata && !params.checkm_metadata && params.gtdb_metadata) {
        ch_gtdbtk_metadata
        .map{ accno, gtdbtk -> [ accno, gtdbtk ] }
        .map { accno, gtdbtk ->
                [
                    accno,
                    [
                    accno: accno[0],
                    checkm_completeness: "",
                    checkm_contamination: "",
                    checkm_strain_heterogeneity: "",
                    contig_count: "",
                    genome_size: "",
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                    ]
                ]
            }
        .set { ch_gtdbtk_metadata }

    ch_gtdbtk_metadata
            .map {
                accno, gtdbtk -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 !in it }
            .map{ it[0] }
            .join(
                ch_gtdbtk_metadata.map { accno, gtdbtk -> [ gtdbtk.accno, gtdbtk ] }
            )
            .map{ it[1]}
        .set{ ch_gtdbtk_filtered }

        ch_gtdbtk_metadata
            .map {
                accno, gtdbtk -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 in it }
            .map{ it[0] }
            .join(
                ch_gtdb_metadata.map { gtdb -> [ gtdb.accno, gtdb ] }
            )
            .map{ it[1]}
            .mix(ch_gtdbtk_filtered)
        .set{ ch_metadata }
    } else if( params.gtdbtk_metadata && params.checkm_metadata && !params.gtdb_metadata) {
        ch_gtdbtk_metadata
            .map{ accno, gtdbtk -> [ accno, gtdbtk ] }
            .join( ch_checkm_metadata
            .map{ accno, checkm -> [ accno, checkm ] }
            )
            .map { accno, gtdbtk, checkm ->
                    [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                    ]
            }
        .set { ch_metadata }
    } else if( !params.gtdbtk_metadata && !params.checkm_metadata && !params.gtdb_metadata) {
        ch_metadata = Channel.empty()
    } else if( !params.gtdbtk_metadata && params.checkm_metadata && !params.gtdb_metadata) {
        ch_checkm_metadata
            .map{ accno, checkm ->
                [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: "",
                    gtdb_representative: "",
                    gtdb_taxonomy: ""
                ]
            }
            .set { ch_metadata }
    } else if( params.gtdbtk_metadata && !params.checkm_metadata && !params.gtdb_metadata) {
        ch_gtdbtk_metadata
        .map{ accno, gtdbtk -> [ accno, gtdbtk ] }
        .map { accno, gtdbtk ->
                [
                    accno: accno[0],
                    checkm_completeness: "",
                    checkm_contamination: "",
                    checkm_strain_heterogeneity: "",
                    contig_count: "",
                    genome_size: "",
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                ]
            }
        .set { ch_metadata }
    } else if( !params.gtdbtk_metadata && !params.checkm_metadata && params.gtdb_metadata) {
        ch_metadata = ch_gtdb_metadata
    }

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
    ch_short_reads = Channel.empty()
    ch_short_reads = CAT_FASTQ.out.reads.mix(ch_short_reads_catskipped)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_short_reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    ch_collect_stats = ch_short_reads.collect { meta, fasta -> meta.id }.map { [ [ id:"magmap" ], it ] }
    if ( params.skip_trimming ) {
        ch_collect_stats
            .map { meta, samples -> [ meta, samples, [] ] }
            .set { ch_collect_stats }

    } else {
        if ( params.se_reads ) {
            ch_collect_stats
                .combine(FASTQC_TRIMGALORE.out.trim_log.collect { meta, report -> report }.map { [ it ] })
                .set { ch_collect_stats }
        } else {
            ch_collect_stats
                .combine(FASTQC_TRIMGALORE.out.trim_log.collect { meta, report -> report[0] }.map { [ it ] })
                .set { ch_collect_stats }
        }
    }

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, params.sequence_filter )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { it[1] }.map { [ it ] }
        ch_versions   = ch_versions.mix(BBMAP_BBDUK.out.versions)
        ch_collect_stats
            .combine(ch_bbduk_logs)
            .set {ch_collect_stats}
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = Channel.empty()
        ch_collect_stats
            .map { [ it[0], it[1], it[2], [] ] }
            .set { ch_collect_stats }
    }

    //
    // SUBWORKFLOW: Use SOURMASH on samples reads and genomes to reduce the number of the latter
    //
    // we create a channel for ncbi genomes only when sourmash is called
    if ( params.sourmash ) {
        SOURMASH(ch_clean_reads, ch_indexes, ch_genomeinfo, ch_genome_infos, params.ksize)
        ch_versions = ch_versions.mix(SOURMASH.out.versions)
        ch_genomes = SOURMASH.out.filtered_genomes
    } else {
        ch_genomeinfo
            .map { [
                accno: it.accno,
                genome_fna: file(it.genome_fna),
                genome_gff: it.genome_gff ? file(it.genome_gff) : ''
                ] }
            .set{ ch_genomes }
    }

    // filter the genomes for the metadata and save it in results/summary_tables directory
    if( params.gtdbtk_metadata || params.checkm_metadata || params.gtdb_metadata) {
        ch_header = Channel
            .of( "accno\tcheckm_completeness\tcheckm_contamination\t \
            checkm_strain_heterogeneity\tcontig_count\tgenome_size\t \
            gtdb_genome_representative\tgtdb_representative\tgtdb_taxonomy")

        ch_metadata
            .map { [ it.accno, it ] }
            .join( ch_genomes
                .map{it.accno}
            )
            .map{ it[1] }
            .map {
                "$it.accno\t$it.checkm_completeness\t \
                $it.checkm_contamination\t$it.checkm_strain_heterogeneity\t \
                $it.contig_count\t$it.genome_size\t \
                $it.gtdb_genome_representative\t$it.gtdb_representative\t \
                $it.gtdb_taxonomy" }
        .set { ch_metadata }

        ch_header
            .concat( ch_metadata )
            .collectFile(name: "magmap.summary_table.taxonomy.tsv",
            newLine: true,
            storeDir: "${params.outdir}/summary_tables")
    }

    //
    // MODULE: Prokka - get gff for all genomes that lack of it
    //
    ch_genomes
        .filter{ !it.genome_gff }
        .map{ [ [id: it.accno ] , it.genome_fna ] }
        .set { ch_no_gff }

    // GUNZIP gff files provided by the user
    ch_genomes
        .filter{ it.genome_gff }
        .set{ ch_genomes_with_gff }

    GUNZIP(ch_no_gff)

    PROKKA(GUNZIP.out.gunzip, [], [])

    ch_genomes_with_gff
        .mix(PROKKA.out.gff
            .map{ meta, gff -> [ meta.id  , [ meta.id, gff ] ] }
            .join(ch_no_gff.map { meta, fna -> [ meta.id , [ meta.id, fna ] ] } )
            .map{ meta, gff, fna -> [ accno: gff[0], genome_fna: fna[1], genome_gff: gff[1] ] })
        .set{ ch_ready_genomes }

    //
    // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
    //
    CREATE_BBMAP_INDEX ( ch_ready_genomes.map{ it.genome_fna } )
    ch_versions = ch_versions.mix(CREATE_BBMAP_INDEX.out.versions)

    //
    // SUBWORKFLOW: Concatenate gff files
    //
    CAT_GFFS ( ch_ready_genomes.map{ it.genome_gff } )
    ch_versions = ch_versions.mix(CAT_GFFS.out.versions)

    //
    // BBMAP ALIGN. Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, CREATE_BBMAP_INDEX.out.index )
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    //
    // SUBWORKFLOW: sort bam file and produce statistics
    //
    BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, CREATE_BBMAP_INDEX.out.genomes_fnas)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(CAT_GFFS.out.gff.map { it[1] })
        .set { ch_featurecounts }

    ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it[1]}.map { [ it ] })
        .set { ch_collect_stats }

    //
    // MODULE: FeatureCounts
    //
    FEATURECOUNTS ( ch_featurecounts )
    ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    FEATURECOUNTS.out.counts
        .collect() { it[1] }
        .map { [ [ id:'all_samples'], it ] }
        .set { ch_collect_features }

    COLLECT_FEATURECOUNTS ( ch_collect_features )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { it[1]}.map { [ it ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { it[1]}
    ch_collect_stats
        .combine(ch_fcs_for_stats)
        .set { ch_collect_stats }

    //
    // Collect statistics from the pipeline
    //
    COLLECT_STATS(ch_collect_stats)
    ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'magmap_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


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

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
