//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../../modules/nf-core/fastqc/'
include { TRIMGALORE       } from '../../../modules/nf-core/trimgalore/'

workflow FASTQC_TRIMGALORE {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc   // boolean: true/false
    skip_trimming // boolean: true/false

    main:
    ch_versions    = Channel.empty()
    fastqc_html    = Channel.empty()
    fastqc_zip     = Channel.empty()

    if (!skip_fastqc) {
        fastqc_html    = FASTQC ( reads ).html
        fastqc_zip     = FASTQC.out.zip
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    trim_reads = reads
    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()

    if (!skip_trimming) {
        trim_reads  = TRIMGALORE ( reads ).reads
        trim_html   = TRIMGALORE.out.html
        trim_zip    = TRIMGALORE.out.zip
        trim_log    = TRIMGALORE.out.log
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]

    versions = ch_versions // channel: [ versions.yml ]

}
