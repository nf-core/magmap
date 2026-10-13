//
// Search proteins with HMM profiles, rank the hits, attach genome accessions and extract hit sequences.
//

include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { HMMER_FORMATTSV } from '../../../modules/nf-core/hmmer/formattsv/main'
include { DUCKDB_TABLE2PARQUET as HMMSEARCH_TABLE2PARQUET } from '../../../modules/nf-core/duckdb/table2parquet/main'
include { HMMER_HMMRANK } from '../../../modules/nf-core/hmmer/hmmrank/main'
include { TIDYVERSE_JOINHMMRANKACCNO } from '../../../modules/local/tidyverse/joinhmmrankaccno/main'
include { SEQTK_HMMHITFAAS } from '../../../modules/local/seqtk/hmmhitfaas/main'

workflow HMMCLASSIFY {
    take:
    ch_hmms // channel: path(hmm)
    ch_faa // channel: [ val(meta), path(faa) ]: proteins to search
    ch_genomes2orfs // channel: path(genomes2orfs): genome accession to ORF index

    main:
    HMMER_HMMSEARCH(
        ch_faa.combine(ch_hmms).map { meta, faa, hmm ->
            def profile = hmm.name.replaceAll(/(\.hmm)?(\.gz)?$/, '')
            [[id: "${meta.id}.${profile}", profile: profile, faa_id: meta.id], hmm, faa, false, true, false]
        }
    )

    // toSortedList emits an empty list when no profiles were given, hence the filter
    HMMER_FORMATTSV(
        HMMER_HMMSEARCH.out.target_summary.map { meta, tbl -> [meta.profile, meta.faa_id, tbl] }.toSortedList { a, b -> a[0] <=> b[0] }.filter { rows -> rows }.map { rows -> [[id: rows[0][1]], rows.collect { r -> r[0] }, rows.collect { r -> r[2] }] },
        'tblout',
    )

    HMMSEARCH_TABLE2PARQUET(HMMER_FORMATTSV.out.tsv)

    HMMER_HMMRANK(HMMSEARCH_TABLE2PARQUET.out.parquet.map { meta, parquet -> [meta, parquet, []] })

    TIDYVERSE_JOINHMMRANKACCNO(HMMER_HMMRANK.out.hmmrank, ch_genomes2orfs)

    SEQTK_HMMHITFAAS(
        HMMER_HMMRANK.out.hmmrank.combine(ch_faa.map { _meta, faa -> faa })
    )

    emit:
    tbl = HMMER_HMMSEARCH.out.target_summary // channel: [ val(meta), path(tbl) ]
    hmmrank = TIDYVERSE_JOINHMMRANKACCNO.out.hmmrank // channel: [ val(meta), path(hmmrank) ]
    faas = SEQTK_HMMHITFAAS.out.faas // channel: [ val(meta), path(faas) ]
}
