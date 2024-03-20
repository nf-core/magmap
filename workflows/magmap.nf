/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMagmap.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local
//
include { COLLECT_FEATURECOUNTS } from '../modules/local/collect_featurecounts'
include { COLLECT_STATS         } from '../modules/local/collect_stats'
include { FILTER_GENOMES        } from '../modules/local/filter_genomes'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { CHECKM_QC   } from '../subworkflows/local/checkm_qc'

//
// SUBWORKFLOW: Local
//
include { FASTQC_TRIMGALORE   } from '../subworkflows/local/fastqc_trimgalore'
include { CAT_GFFS            } from '../subworkflows/local/concatenate_gff'
include { CREATE_BBMAP_INDEX  } from '../subworkflows/local/create_bbmap_index'
include { SOURMASH            } from '../subworkflows/local/sourmash'
include { ARIA2_UNTAR         } from '../subworkflows/local/aria2_untar'
include { GTDBTK              } from '../subworkflows/local/gtdbtk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { BBMAP_BBDUK                            } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_ALIGN                            } from '../modules/nf-core/bbmap/align/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'
include { GUNZIP                                 } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFFS                  } from '../modules/nf-core/gunzip/main'
include { PROKKA                                 } from '../modules/nf-core/prokka/main'
//include { ARIA2                                  } from '../modules/nf-core/aria2/main'
//include { UNTAR                                  } from '../modules/nf-core/untar/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { BAM_SORT_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

gtdb = ( params.skip_binqc || params.skip_gtdbtk ) ? false : params.gtdb_db

if (gtdb) {
    gtdb = file( "${gtdb}", checkIfExists: true)
    gtdb_mash = params.gtdb_mash ? file("${params.gtdb_mash}", checkIfExists: true) : []
} else {
    gtdb = []
}

// Info required for completion email and summary
def multiqc_report = []

workflow MAGMAP {

    ch_versions = Channel.empty()

    if ( !params.skip_binqc ) {
        ARIA2_UNTAR()
        ch_checkm_db = ARIA2_UNTAR.out.checkm_db
        ch_versions  = ARIA2_UNTAR.out.versions
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // INPUT: if user provides, populate ch_genomeinfo with a table that provides the genomes to filter with sourmash
    //
    if ( params.genomeinfo) {
        Channel
            .fromPath( params.genomeinfo )
            .splitCsv( sep: ',', header: true )
            .set { ch_genomeinfo }
    }

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
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    ch_collect_stats = INPUT_CHECK.out.reads.collect { it[0].id }.map { [ [ id: "magmap" ], it ] }
    if ( params.skip_trimming ) {
        ch_collect_stats
            .map { [ it[0], it[1], [] ] }
            .set { ch_collect_stats }
    } else {
        ch_collect_stats
            .combine(FASTQC_TRIMGALORE.out.trim_log.collect { it[1][0] }.map { [ it ] })
            .set { ch_collect_stats }
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
        SOURMASH(ch_clean_reads, ch_indexes, ch_genomeinfo, ch_genome_infos)
        ch_versions = ch_versions.mix(SOURMASH.out.versions)
        ch_genomes = SOURMASH.out.filtered_genomes
    } else {
        ch_genomes = ch_genomeinfo
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
        .map { [ [id: it.accno], it.genome_gff ] }
        .set{ gff_to_gunzip }

     GUNZIP_GFFS(gff_to_gunzip)
     GUNZIP_GFFS.out.gunzip
         .map{ meta, gff -> [ [id: meta.id], gff ] }
         .join(ch_genomes
             .filter{ it.genome_gff }
             .map { [ [id:it.accno], it.genome_fna ] })
         .map{ meta, gff, fna -> [ accno: meta.id, genome_fna: fna, genome_gff: gff ] }
         .set { ch_genomes_gunzipped_gff }

    GUNZIP(ch_no_gff)

    PROKKA(GUNZIP.out.gunzip, [], [])

    ch_genomes_gunzipped_gff
        .mix(PROKKA.out.gff
            .ifEmpty([])
            .map{ meta, gff -> [ meta.id  , [ meta.id, gff ] ] }
            .join(ch_no_gff.map { meta, fna -> [ meta.id , [ meta.id, fna ] ] } )
            .map{ meta, gff, fna -> [ accno: gff[0], genome_fna: fna[1], genome_gff: gff[1] ] })
        .set{ ch_ready_genomes }

    //
    // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
    //
    def i = 0

    ch_ready_genomes
    .map{ it.genome_fna }
    .flatten()
    .collate(1000)
    .map{ [ [ id: "all_references${i++}" ], it ] }
    .set { ch_genomes_fnas }

    CREATE_BBMAP_INDEX ( ch_genomes_fnas )
    ch_versions = ch_versions.mix(CREATE_BBMAP_INDEX.out.versions)

    //
    // CheckM
    //
    if (!params.skip_binqc){
        CHECKM_QC (
           ch_genomes_fnas.groupTuple(),
           ch_checkm_db.map { meta, db -> db }
        )
        ch_checkm_summary = CHECKM_QC.out.summary
       ch_versions       = ch_versions.mix(CHECKM_QC.out.versions)
    }

    //
    // GTDB-tk: taxonomic classifications using GTDB reference
    //
    if ( !params.skip_gtdbtk ) {
        ch_gtdbtk_summary = Channel.empty()
        if ( gtdb ){
            GTDBTK (
            ch_genomes_fnas,
            ch_checkm_summary,
            gtdb,
            gtdb_mash
            )
        ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
        ch_gtdbtk_summary = GTDBTK.out.summary
        }
    } else {
        ch_gtdbtk_summary = Channel.empty()
    }

    //
    // SUBWORKFLOW: Concatenate gff files
    //
    CAT_GFFS ( ch_ready_genomes.map{ [ [id: "gffs"], it.genome_gff ] } )
    ch_versions = ch_versions.mix(CAT_GFFS.out.versions)

    //
    // BBMAP ALIGN. Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, CREATE_BBMAP_INDEX.out.index )
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    // //
    // // SUBWORKFLOW: sort bam file and produce statistics
    // //
    // BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, CREATE_BBMAP_INDEX.out.genomes_fnas )
    // ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // BAM_SORT_STATS_SAMTOOLS.out.bam
    //     .combine(CAT_GFFS.out.gff.map { it[1] })
    //     .set { ch_featurecounts }

    // ch_collect_stats
    //     .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it[1]}.map { [ it ] })
    //     .set { ch_collect_stats }

    // //
    // // MODULE: FeatureCounts
    // //
    // FEATURECOUNTS ( ch_featurecounts )
    // ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)

    // //
    // // MODULE: Collect featurecounts output counts in one table
    // //
    // FEATURECOUNTS.out.counts
    //     .collect() { it[1] }
    //     .map { [ [ id:'all_samples'], it ] }
    //     .set { ch_collect_feature }

    // COLLECT_FEATURECOUNTS ( ch_collect_feature )
    // ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    // ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { it[1]}.map { [ it ] }
    // ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { it[1]}
    // ch_collect_stats
    //     .combine(ch_fcs_for_stats)
    //     .set { ch_collect_stats }

    // //
    // // Collect statistics from the pipeline
    // //
    // COLLECT_STATS(ch_collect_stats)
    // ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    //
    // MODULE: custom dump software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMagmap.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMagmap.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
