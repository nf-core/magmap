include { KOFAMSCAN_DOWNLOAD          } from '../../../modules/local/kofamscan/download'
include { KOFAMSCAN as KOFAMSCAN_SCAN } from '../../../modules/nf-core/kofamscan'
include { KOFAMSCAN_TABLE             } from '../../../modules/local/kofamscan/table'
include { KOFAMSCAN_VERSION           } from '../../../modules/local/kofamscan/version'

workflow KOFAMSCAN {

    take:
    ch_proteins     // channel: [ val(meta), path(faa) ], one per genome; meta.id is the accession
    ko_list_url     //  string: URL or path of the KOfam ko_list.gz file
    profiles_url    //  string: URL or path of the KOfam profiles.tar.gz archive

    main:

    KOFAMSCAN_DOWNLOAD( file(ko_list_url), file(profiles_url) )

    KOFAMSCAN_SCAN(
        ch_proteins,
        KOFAMSCAN_DOWNLOAD.out.koprofiles,
        KOFAMSCAN_DOWNLOAD.out.ko_list
    )

    KOFAMSCAN_TABLE(
        KOFAMSCAN_SCAN.out.tsv
            .map { _meta, tsv -> tsv }
            .collect()
            .map { tsvs -> [ [ id: 'magmap' ], tsvs ] }
    )

    KOFAMSCAN_VERSION()

    emit:
    tsv = KOFAMSCAN_TABLE.out.tsv // channel: [ val(meta), path(tsv.gz) ]
}
