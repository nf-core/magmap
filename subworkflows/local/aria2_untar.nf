/*
 * CheckM: Quantitative measures for the assessment of genome assembly
 */

include { ARIA2                             } from '../../modules/nf-core/aria2/main'
include { UNTAR                             } from '../../modules/nf-core/untar/main'

workflow ARIA2_UNTAR {
    take:

    main:
    ch_versions = Channel.empty()

    ARIA2 ([ [id:"download_checkm_db"] , params.checkm_download_url])
        ch_versions = ch_versions.mix(ARIA2.out.versions)
        UNTAR(ARIA2.out.downloaded_file)
        ch_checkm_db =
        ch_versions = ch_versions.mix(UNTAR.out.versions)

    emit:
    checkm_db = UNTAR.out.untar
    versions     = ch_versions
}
