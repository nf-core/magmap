//
// Subworkflow with functionality specific to the nf-core/magmap pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version                 // boolean: Display version and exit
    validate_params         // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs        // boolean: Do not use coloured log outputs
    nextflow_cli_args       //   array: List of positional nextflow CLI args
    outdir                  //  string: The output directory where the results will be saved
    input                   //  string: Path to input samplesheet
    help                    // boolean: Display help message and exit
    help_full               // boolean: Show the full help message
    show_hidden             // boolean: Show hidden parameters in the help message
    genomeinfo              //  string: Path to user-provided genome sheet
    remote_genome_sources   //  string: Comma-separated list of NCBI-style genome summary files
    genome_store_dir        //  string: Path to a directory where genome annotation files will be stored
    prokka_store_dir        //  string: Path to a directory where Prokka annotation files will be stored
    annotator               //  string: 'prokka', 'bakta_supported_only' or 'bakta_all' -- which tool(s) to annotate genomes lacking a GFF with
    bakta_db                //  string: Path to a directory where the Bakta database is (or will be) stored
    bakta_store_dir         //  string: Path to a directory where Bakta annotation files will be stored
    indexes                 //  string: Path to user-provided Sourmash index file
    genomeset_mode          //  string: Genomeset mode ('sample' or 'joint')
    species_preference      //  string: 'all', 'local', 'completeness' or 'gtdb' to indicate preferred genome for a species
    gtdb_metadata           //  string: Paths to GTDB metadata files
    gtdbtk_metadata         //  string: Path to GTDB-Tk metadata file
    checkm_metadata         //  string: Path to GTDB metadata file
    features                //  string: Comma-separated string of feature types

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/magmap ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/magmap/blob/master/CITATIONS.md
"""
    if (_monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create a channel from input file provided through input
    //

    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    //
    // INPUT: if the user provides --genomeinfo, populate ch_genomeinfo with a table that provides the genomes to filter with sourmash
    //
    ch_genomeinfo = channel.empty()
    if ( genomeinfo ) {
        channel
            .fromPath(genomeinfo)
            .splitCsv(sep: ',', header: true)
            .map { it -> [
                    accno: it.accno,
                    genome_fna: file(it.genome_fna),
                    genome_gff: it.genome_gff ? file(it.genome_gff) : []
                ]
            }
            .set { ch_genomeinfo }
    }

    //
    // INPUT: genome info from ncbi
    //
    ch_remote_genome_sources = channel.empty()
    if ( remote_genome_sources ) {
        ch_remote_genome_sources = channel
            .of(remote_genome_sources.split(','))
            .map { it -> file(it) }
    }

    //
    // Make sure that the directories for genome and annotation storage exists
    //
    if ( genome_store_dir ) {
        d = new File("${genome_store_dir}")
        if ( ! d.exists() ) { d.mkdirs() }
    }
    if ( prokka_store_dir ) {
        d = new File("${prokka_store_dir}")
        if ( ! d.exists() ) { d.mkdirs() }
    }
    if ( bakta_db && annotator != 'prokka' ) {
        d = new File("${bakta_db}")
        if ( ! d.exists() ) { d.mkdirs() }
    }
    if ( bakta_store_dir && annotator != 'prokka' ) {
        d = new File("${bakta_store_dir}")
        if ( ! d.exists() ) { d.mkdirs() }
    }

    //
    // INPUT: if the user provides, populate ch_indexes
    //
    ch_indexes = channel.empty()
    if ( indexes ) {
        ch_indexes = channel.fromPath(indexes.tokenize(','))
    }

    //
    // Return error if user asks for sourmash filtering but doesn't provide indexes.
    //
    if (genomeset_mode == 'sample' && !indexes) {
        error("You have asked to run sourmash sample filtering but have not provided any Sourmash indexes. Please provide --indexes or set --skip_sourmash to true.")
    }

    //
    // Return an error if the user asks for a species_preference that requires genome
    // metadata we don't have.
    //
    if (species_preference in ['local', 'completeness', 'gtdb'] && (!gtdb_metadata || !gtdbtk_metadata)) {
        error("--species_preference '${species_preference}' requires genome metadata. Please provide both --gtdb_metadata and --gtdbtk_metadata, or set --species_preference to 'all'.")
    }
    if (species_preference in ['completeness', 'gtdb'] && !checkm_metadata) {
        error("--species_preference '${species_preference}' additionally requires --checkm_metadata.")
    }

    //
    // Take care of genome metadata files
    //
    ch_gtdb_metadata = channel.empty()
    if ( gtdb_metadata ) {
        ch_gtdb_metadata = channel
            .of(gtdb_metadata.split(','))
            .map { it -> file(it) }
    }

    ch_gtdbtk_metadata = channel.empty()
    if ( gtdbtk_metadata ) {
        ch_gtdbtk_metadata = channel
            .of(gtdbtk_metadata.split(','))
            .map { it -> file(it) }
    }

    ch_checkm_metadata = channel.empty()
    if ( checkm_metadata ) {
        ch_checkm_metadata = channel
            .of(checkm_metadata.split(','))
            .map { it -> file(it) }
    }

    ch_features = channel.of(
        ['CDS'] + features.split(','))
        .flatten()
        .unique()

    emit:
    samplesheet             = ch_samplesheet
    genomeinfo              = ch_genomeinfo
    remote_genome_sources   = ch_remote_genome_sources
    indexes                 = ch_indexes
    gtdb_metadata           = ch_gtdb_metadata
    gtdbtk_metadata         = ch_gtdbtk_metadata
    checkm_metadata         = ch_checkm_metadata
    features                = ch_features
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
// Not implemented, since we don't have a genome param
//
def getGenomeAttribute(_attribute) {
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
// Not implemented, since we don't have a genome param
//
def genomeExistsError() {
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Kept in sync with CITATIONS.md's "Pipeline tools" section.
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "Trim Galore!,",
            "sourmash (Brown & Irber 2016),",
            "Prokka (Seemann 2014),",
            "Bakta (Schwengers et al. 2021),",
            "gffread (Pertea & Pertea 2020),",
            "BBMap,",
            "Samtools (Li et al. 2009),",
            "GTDB-Tk (Chaumeil et al. 2020),",
            "CheckM (Parks et al. 2015),",
            "featureCounts (Liao et al. 2014),",
            "R (R Core Team 2025),",
            "Tidyverse (Wickham et al. 2019),",
            "data.table (Barrett et al. 2025),",
            "DuckDB (Raasveldt & Mühleisen 2019)",
            "and MultiQC (Ewels et al. 2016)."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Kept in sync with CITATIONS.md's "Pipeline tools" section.
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Trim Galore!, URL: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/.</li>",
            "<li>Brown, C.T., Irber Junior, L. C. Sourmash: a library for MinHash sketching of DNA. 2016, The Journal of Open Source Software. DOI: 10.21105/joss.00027</li>",
            "<li>Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063</li>",
            "<li>Schwengers O., Jelonek L., Dieckmann M. A., Beyvers S., Blom J., Goesmann A. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 2021;7(11):000685. doi: 10.1099/mgen.0.000685. PMID: 34739369; PMCID: PMC8743544.</li>",
            "<li>Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Research 2020, 9:304. doi: 10.12688/f1000research.23297.2.</li>",
            "<li>BBMap, URL: https://sourceforge.net/projects/bbmap/.</li>",
            "<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PMID: 19505943; PMCID: PMC2723002.</li>",
            "<li>Pierre-Alain Chaumeil, Aaron J Mussig, Philip Hugenholtz, Donovan H Parks, GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database, Bioinformatics, Volume 36, Issue 6, March 2020, Pages 1925–1927</li>",
            "<li>Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Res. 2015 Jul;25(7):1043-55. doi: 10.1101/gr.186072.114. Epub 2015 May 14. PMID: 25977477; PMCID: PMC4484387.</li>",
            "<li>Liao Y, Smyth GK and Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47(8):e47, 2019.</li>",
            "<li>Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014.</li>",
            "<li>Liao Y, Smyth GK and Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108, 2013.</li>",
            "<li>R Core Team (2025): R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.</li>",
            "<li>Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019): Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686</li>",
            "<li>Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T, Schwendinger B, Krylov I (2025): data.table: Extension of `data.frame`. doi:10.32614/CRAN.package.data.table</li>",
            "<li>Raasveldt M, Mühleisen H. DuckDB: an Embeddable Analytical Database. In: Proceedings of the 2019 International Conference on Management of Data (SIGMOD '19). 2019 Jun 25:1981-1984. doi: 10.1145/3299869.3320212.</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
