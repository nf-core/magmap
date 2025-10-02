#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/magmap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/magmap
    Website: https://nf-co.re/magmap
    Slack  : https://nfcore.slack.com/channels/magmap
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAGMAP                  } from './workflows/magmap'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_magmap_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_magmap_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_MAGMAP {

    take:
    samplesheet                 // channel: samplesheet read in from --input
    genomeinfo                  // channel: genome information sheet read in from --genomeinfo
    remote_genome_sources       // channel: NCBI-style genome summary files read in via --remote_genome_sources
    indexes                     // channel: user-provided Sourmash indexes
    sequence_filter             // channel: fasta file for BBDuk
    gtdb_metadata               // channel: GTDB metadata files
    gtdbtk_metadata             // channel: GTDB-Tk metadata files
    checkm_metadata             // channel: CheckM/CheckM2 metadata files
    skip_kraken2                // boolean: run Kraken2 or not
    kraken2_db                  // string: path to Kraken2 database
    kraken2_db_url              // string: URL to download Kraken2 database
    skip_sourmash               // boolean: skip Sourmash or not
    sourmash_ksize              // integer
    sourmash_save_unassigned    // boolean
    sourmash_save_matches_sig   // boolean
    sourmash_save_prefetch      // boolean
    sourmash_save_prefetch_csv  // boolean
    features                    // channel: types of features to call
    skip_fastqc                 // boolean
    skip_qc                     // boolean
    skip_trimming               // boolean
    outdir                      //  string: path to output directory

    main:

    //
    // WORKFLOW: Run pipeline
    //
    MAGMAP (
        samplesheet,
        genomeinfo,
        remote_genome_sources,
        indexes,
        sequence_filter,
        gtdb_metadata,
        gtdbtk_metadata,
        checkm_metadata,
        skip_kraken2,
        kraken2_db,
        kraken2_db_url,
        skip_sourmash,
        sourmash_ksize,
        sourmash_save_unassigned,
        sourmash_save_matches_sig,
        sourmash_save_prefetch,
        sourmash_save_prefetch_csv,
        features,
        skip_fastqc,
        skip_qc,
        skip_trimming,
        outdir
    )
    emit:
    multiqc_report = MAGMAP.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.genomeinfo,
        params.remote_genome_sources,
        params.kraken2_store_dir,
        params.genome_store_dir,
        params.indexes,
        params.gtdb_metadata,
        params.gtdbtk_metadata,
        params.checkm_metadata,
        params.features
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_MAGMAP (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.genomeinfo,
        PIPELINE_INITIALISATION.out.remote_genome_sources,
        PIPELINE_INITIALISATION.out.indexes,
        params.sequence_filter,
        PIPELINE_INITIALISATION.out.gtdb_metadata,
        PIPELINE_INITIALISATION.out.gtdbtk_metadata,
        PIPELINE_INITIALISATION.out.checkm_metadata,
        params.skip_kraken2,
        params.kraken2_db,
        params.kraken2_db_url,
        params.skip_sourmash,
        params.sourmash_ksize,
        params.sourmash_save_unassigned,
        params.sourmash_save_matches_sig,
        params.sourmash_save_prefetch,
        params.sourmash_save_prefetch_csv,
        PIPELINE_INITIALISATION.out.features,
        params.skip_fastqc,
        params.skip_qc,
        params.skip_trimming,
        params.outdir
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_MAGMAP.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
