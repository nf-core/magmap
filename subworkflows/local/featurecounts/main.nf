//
// Create a file that summarises the feature counts for each feature type (CDS, rRNA, tRNA, tmRNA) and the versions of the software used.
//

include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS   } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_RRNA  } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_TRNA  } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_TMRNA } from '../../../modules/nf-core/subread/featurecounts/main'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_CDS   } from '../../../modules/local/collect_featurecounts'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_RRNA  } from '../../../modules/local/collect_featurecounts'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_TRNA  } from '../../../modules/local/collect_featurecounts'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_TMRNA } from '../../../modules/local/collect_featurecounts'
    
 workflow FEATURECOUNTS {   

    take:
        ch_featurecounts

    main:
    ch_versions = Channel.empty()
    ch_cds_counts = Channel.empty()

    if ( params.count_cds ) {
        FEATURECOUNTS_CDS ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_CDS.out.versions)

        FEATURECOUNTS_CDS.out.counts
        .collect() { it[1] }
        .map { [ [ id:'all_samples_cds'], it ] }
        .set { ch_collect_features_cds }

        COLLECT_FEATURECOUNTS_CDS   ( ch_collect_features_cds )
        ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTS_CDS.out.versions)
        ch_cds_counts = COLLECT_FEATURECOUNTS_CDS.out.counts
    }

    ch_rrna_counts = Channel.empty()
    if ( params.count_rrna ) {
        FEATURECOUNTS_RRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_RRNA.out.versions)

        FEATURECOUNTS_RRNA.out.counts
        .collect() { it[1] }
        .map { [ [ id:'all_samples_rrna'], it ] }
        .set { ch_collect_features_rrna }

        COLLECT_FEATURECOUNTS_RRNA  ( ch_collect_features_rrna)
        ch_rrna_counts = COLLECT_FEATURECOUNTS_RRNA.out.counts
    }

    ch_trna_counts = Channel.empty()
    if ( params.count_trna ) {
        FEATURECOUNTS_TRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_TRNA.out.versions)

        FEATURECOUNTS_TRNA.out.counts
        .collect() { it[1] }
        .map { [ [ id:'all_samples_trna'], it ] }
        .set { ch_collect_features_trna }

        COLLECT_FEATURECOUNTS_TRNA  ( ch_collect_features_trna )
        ch_trna_counts = COLLECT_FEATURECOUNTS_TRNA.out.counts
    }

    ch_tmrna_counts = Channel.empty()
    if ( params.count_tmrna ) {
        FEATURECOUNTS_TMRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_TMRNA.out.versions)

        FEATURECOUNTS_TMRNA.out.counts
        .collect() { it[1] }
        .map { [ [ id:'all_samples_tmrna'], it ] }
        .set { ch_collect_features_tmrna }

        COLLECT_FEATURECOUNTS_TMRNA ( ch_collect_features_tmrna )
        ch_tmrna_counts = COLLECT_FEATURECOUNTS_TMRNA.out.counts
    }

    // Collect feature counts into one channel for summary
    ch_fcs = Channel.empty()
    ch_fcs = ch_fcs.mix(
        ch_cds_counts,
        ch_rrna_counts,
        ch_trna_counts,
        ch_tmrna_counts
    ).collect { it[1]}.map { [ it ] }

    emit:
        collected_features = ch_fcs
        versions           = ch_versions

 }