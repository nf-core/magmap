/*
 * CheckM: Quantitative measures for the assessment of genome assembly
 */

include { CHECKM_QA                         } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF                  } from '../../modules/nf-core/checkm/lineagewf/main'
include { COMBINE_TSV as COMBINE_CHECKM_TSV } from '../../modules/local/combine_tsv'

workflow CHECKM_QC {
    take:
    bins       // channel: [ val(meta), path(bin) ]
    checkm_db

    main:
    ch_versions = Channel.empty()

    ch_input_checkmdb = checkm_db ? checkm_db : []

    CHECKM_LINEAGEWF ( bins, "fna.gz", checkm_db )
    ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

    ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
        .join(CHECKM_LINEAGEWF.out.marker_file)
        .map{
            meta, dir, marker ->
            [ meta, dir, marker, []]
        }

    CHECKM_QA ( ch_checkmqa_input, [] )
    ch_versions = ch_versions.mix(CHECKM_QA.out.versions.first())

    COMBINE_CHECKM_TSV ( CHECKM_QA.out.output.map{it[1]}.collect() )
    ch_versions = ch_versions.mix(COMBINE_CHECKM_TSV.out.versions.first())

    emit:
    summary    = COMBINE_CHECKM_TSV.out.combined
    checkm_tsv = CHECKM_QA.out.output
    versions   = ch_versions
}
