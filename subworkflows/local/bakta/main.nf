//
// Annotate genomes with Bakta: download (or reuse a cached/user-supplied) database, then
// run Bakta on genomes lacking a GFF. Bakta's own version reporting is handled by the
// storeDir-decoupled BAKTA_VERSION process -- see modules/local/bakta_version.
//

include { BAKTA_BAKTADBDOWNLOAD } from '../../../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA           } from '../../../modules/nf-core/bakta/bakta/main'
include { BAKTA_VERSION         } from '../../../modules/local/bakta_version/main'

workflow BAKTA {

    take:
    ch_fasta // channel: [ val(meta), path(fasta) ] -- genomes without a GFF to annotate with Bakta

    main:
    // Only download the (potentially large) Bakta database if there's actually a genome
    // to annotate with it -- ch_fasta can be empty even when --annotator requests Bakta,
    // e.g. if domain classification routed every genome lacking a GFF to Prokka instead.
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
    gff = BAKTA_BAKTA.out.gff // channel: [ val(meta), path(gff3) ]
    txt = BAKTA_BAKTA.out.txt // channel: [ val(meta), path(txt) ]
}
