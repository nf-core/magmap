//
// Annotate genomes lacking a GFF with Bakta, downloading or reusing a database.
// BAKTA_VERSION reports the Bakta version, decoupled from the storeDir-cached process.
//

include { BAKTA_BAKTADBDOWNLOAD } from '../../../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA           } from '../../../modules/nf-core/bakta/bakta/main'
include { BAKTA_VERSION         } from '../../../modules/local/bakta_version/main'

workflow BAKTA {

    take:
    ch_fasta // channel: [ val(meta), path(fasta) ]: genomes without a GFF to annotate with Bakta

    main:
    // Download the large database only if there is a genome to annotate; ch_fasta can be
    // empty when --annotator requests Bakta but every genome was routed to Prokka.
    BAKTA_BAKTADBDOWNLOAD(ch_fasta.count().filter { n -> n > 0 })

    BAKTA_BAKTA(
        ch_fasta,
        BAKTA_BAKTADBDOWNLOAD.out.db,
        [],
        [],
        [],
        []
    )

    BAKTA_VERSION()

    emit:
    fna = BAKTA_BAKTA.out.fna // channel: [ val(meta), path(fna) ]
    faa = BAKTA_BAKTA.out.faa // channel: [ val(meta), path(faa) ]
    gff = BAKTA_BAKTA.out.gff // channel: [ val(meta), path(gff3) ]
    txt = BAKTA_BAKTA.out.txt // channel: [ val(meta), path(txt) ]
}
