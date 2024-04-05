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
include { validateInputSamplesheet       } from '../subworkflows/local/utils_nfcore_magmap_pipeline'

//
// SUBWORKFLOW: Local
//
include { FASTQC_TRIMGALORE       } from '../subworkflows/local/fastqc_trimgalore'
include { CAT_GFFS                } from '../subworkflows/local/concatenate_gff'
include { CREATE_BBMAP_INDEX      } from '../subworkflows/local/create_bbmap_index'
include { SOURMASH                } from '../subworkflows/local/sourmash'
include { ARIA2_UNTAR             } from '../subworkflows/local/aria2_untar'
include { GTDBTK                  } from '../subworkflows/local/gtdbtk'
include { PIPELINE_INITIALISATION } from '../subworkflows/local/utils_nfcore_magmap_pipeline'
include { PIPELINE_COMPLETION     } from '../subworkflows/local/utils_nfcore_magmap_pipeline'

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
include { BBMAP_BBDUK                            } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_ALIGN                            } from '../modules/nf-core/bbmap/align/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'
include { GUNZIP                                 } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFFS                  } from '../modules/nf-core/gunzip/main'
include { PROKKA                                 } from '../modules/nf-core/prokka/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { paramsSummaryMap                       } from 'plugin/nf-validation'
include { fromSamplesheet                        } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_SORT_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { UTILS_NEXTFLOW_PIPELINE                } from '../subworkflows/nf-core/utils_nextflow_pipeline/main'
include { UTILS_NFCORE_PIPELINE                  } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { UTILS_NFVALIDATION_PLUGIN              } from '../subworkflows/nf-core/utils_nfvalidation_plugin/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// set an empty multiqc channel
ch_multiqc_files = Channel.empty()

gtdb = ( params.skip_binqc || params.skip_gtdbtk ) ? false : params.gtdb_db

if (gtdb) {
    gtdb = file( "${gtdb}", checkIfExists: true)
    gtdb_mash = params.gtdb_mash ? file("${params.gtdb_mash}", checkIfExists: true) : []
} else {
    gtdb = []
}

workflow MAGMAP {

    take:
    ch_samplesheet       // channel: path(sample_sheet.csv)
    ch_versions          // channel: [ path(versions.yml) ]

    main:
    ch_versions = Channel.empty()

    if ( !params.skip_binqc ) {
        ARIA2_UNTAR()
        ch_checkm_db = ARIA2_UNTAR.out.checkm_db
        ch_versions  = ARIA2_UNTAR.out.versions
    }

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
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .set { ch_fastq }

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    ch_collect_stats = ch_fastq.collect { it[0].id }.map { [ [ id: "magmap" ], it ] }
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
        ch_genomeinfo
            .map { [
                accno: it.accno,
                genome_fna: file(it.genome_fna),
                genome_gff: file(it.genome_gff)
                ] }
            .set{ ch_genomes }
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
        .map { [ [id: it.accno], file(it.genome_gff) ] }
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
        //.mix(PROKKA.out.gff
        //    .ifEmpty([])
        //    .map{ meta, gff -> [ meta.id  , [ meta.id, gff ] ] }
        //    .join(ch_no_gff.map { meta, fna -> [ meta.id , [ meta.id, fna ] ] } )
        //    .map{ meta, gff, fna -> [ accno: gff[0], genome_fna: fna[1], genome_gff: gff[1] ] })
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

    //
    // SUBWORKFLOW: sort bam file and produce statistics
    //
    BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, CREATE_BBMAP_INDEX.out.genomes_fnas )
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
        .set { ch_collect_feature }

    COLLECT_FEATURECOUNTS ( ch_collect_feature )
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

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_magmap_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()
    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{ meta, zip -> zip })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{ meta, idxstats -> idxstats })
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.collect{ meta, summary -> summary })


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_multiqc_report = MULTIQC.out.report.toList()

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
