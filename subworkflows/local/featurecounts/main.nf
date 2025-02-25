//
// Create a file that summarises the feature counts for each feature type (CDS, rRNA, tRNA, tmRNA) and the versions of the software used.
//

include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS   } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_RRNA  } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_TRNA  } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_TMRNA } from '../../../modules/nf-core/subread/featurecounts/main'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_CDS   } from '../../modules/local/collect_featurecounts'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_RRNA  } from '../../modules/local/collect_featurecounts'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_TRNA  } from '../../modules/local/collect_featurecounts'
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_TMRNA } from '../../modules/local/collect_featurecounts'
    
 workflow FEATURECOUNTS {   

    take:
        ch_samples_reads
        ch_indexes
        ch_user_genomeinfo
        ch_ncbi_genomeinfo_files
        ksize

    main:
    ch_cds_counts = Channel.empty()
    if ( params.count_cds ) {
        FEATURECOUNTS_CDS ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_CDS.out.versions)

        COLLECT_FEATURECOUNTS_CDS   ( FEATURECOUNTS_CDS.out.counts.collect   { it[1] } )
        ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTS_CDS.out.versions)
        ch_cds_counts = COLLECT_FEATURECOUNTS_CDS.out.counts
    }

    ch_rrna_counts = Channel.empty()
    if ( params.count_rrna ) {
        FEATURECOUNTS_RRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_RRNA.out.versions)

        COLLECT_FEATURECOUNTS_RRNA  ( FEATURECOUNTS_RRNA.out.counts.collect  { it[1] } )
        ch_rrna_counts = COLLECT_FEATURECOUNTS_RRNA.out.counts
    }

    ch_trna_counts = Channel.empty()
    if ( params.count_trna ) {
        FEATURECOUNTS_TRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_TRNA.out.versions)

        COLLECT_FEATURECOUNTS_TRNA  ( FEATURECOUNTS_TRNA.out.counts.collect  { it[1] } )
        ch_trna_counts = COLLECT_FEATURECOUNTS_TRNA.out.counts
    }

    ch_tmrna_counts = Channel.empty()
    if ( params.count_tmrna ) {
        FEATURECOUNTS_TMRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_TMRNA.out.versions)

        COLLECT_FEATURECOUNTS_TMRNA ( FEATURECOUNTS_TMRNA.out.counts.collect { it[1] } )
        ch_tmrna_counts = COLLECT_FEATURECOUNTS_TMRNA.out.counts
    }

    // Collect feature counts into one channel for summary
    ch_fcs = Channel.empty()
    ch_fcs = ch_fcs.mix(
        ch_cds_counts,
        ch_rrna_counts,
        ch_trna_counts,
        ch_tmrna_counts
    ).collect()

    emit:
        collected_features = ch_fcs
        versions         = ch_versions

 }