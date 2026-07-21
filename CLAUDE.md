# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

nf-core/magmap is a Nextflow (DSL2) pipeline that maps short reads to (potentially large)
collections of prokaryotic genomes and quantifies features (genes) per genome. It's aimed
at metagenomes/metatranscriptomes, as an alternative to single-reference pipelines like
nf-core/rnaseq. Core steps: FastQC -> Trim Galore! -> BBDuk (optional contaminant filter) ->
Sourmash (optional genome selection) -> Prokka/GFF annotation -> BBMap index+align ->
featureCounts -> R-based summary tables -> MultiQC.

## Commands

Requires Nextflow, nf-core tools, and nf-test.

```bash
# Run the pipeline against test data (fast, uses conf/test.config)
nextflow run . -profile test,docker --outdir <OUTDIR>

# Run the full nf-test suite (mirrors CI)
nf-test test --tag test --profile +docker --verbose

# Run a single nf-test file
nf-test test tests/sourmash.nf.test --profile +docker --verbose

# Update snapshots after intentional output changes
nf-test test --tag test --profile +docker --verbose --update-snapshots

# Lint (nf-core conventions)
nf-core pipelines lint .

# Regenerate nextflow_schema.json after changing params in nextflow.config
nf-core pipelines schema build
```

Other test profiles live in `conf/test_*.config` and each has a matching `.nf.test` /
`.nf.test.snap` pair in `tests/` (e.g. `test_sourmash_mix_dupl_species.config` ↔
`tests/sourmash_mix_dupl_species.nf.test`). When adding a new test scenario, add both.

`.nftignore` in `tests/` lists output paths excluded from content-snapshot comparison
(non-deterministic files).

## Architecture

Entry point `main.nf` wires together, in order:
`PIPELINE_INITIALISATION` (subworkflows/local/utils_nfcore_magmap_pipeline) → the
`NFCORE_MAGMAP` wrapper workflow → `MAGMAP` (workflows/magmap.nf, the actual pipeline
logic) → `PIPELINE_COMPLETION`. Parameter validation, samplesheet parsing (via
`plugin/nf-schema`), and genome/index channel setup all happen in
`PIPELINE_INITIALISATION`; look there first when tracing how a CLI param becomes a channel.

### Genome selection is the core complexity

Unlike a typical single-reference nf-core pipeline, magmap builds its set of reference
genomes dynamically from up to three sources, then optionally filters them per-sample:

1. **User-provided genomes** (`--genomeinfo` CSV: `accno,genome_fna,genome_gff`).
2. **Local Sourmash indexes** built from those user genomes (skipped if `--skip_sourmash`).
3. **Remote genomes** looked up in NCBI assembly summary files (`--remote_genome_sources`,
   defaults to RefSeq+GenBank summaries) and pulled in via user-supplied Sourmash indexes
   (`--indexes`), then fetched with `WGET`.

This logic lives in `subworkflows/local/sourmash/main.nf` (workflow `SOURMASH`). It emits
two channels: `joint_filtered_genomes` (one genome set for the whole run) and
`sample_filtered_genomes` (per-sample sets). Which one is actually used downstream depends
on `--genomeset_mode` (`joint` or `sample`), branched in `workflows/magmap.nf` — `joint`
builds one shared BBMap index via `CREATE_BBMAP_INDEX`; `sample` builds a distinct index
per sample from `sample_filtered_genomes` and joins reads to their matching index before
`BBMAP_ALIGN`. `--species_preference` (`all`/`local`/`completeness`/`gtdb`) controls how
duplicate-species genomes are resolved and is consumed later by
`modules/local/tidyverse/joinmetadata`.

Genomes lacking a GFF (mostly ones fetched from NCBI) are annotated with Prokka
(`PROKKA` + `PROKKAGFF2TSV` + `CATPROKKATSVS`); genomes that already have a GFF skip
straight to `PROKKAGFF2TSV`. Both paths reconverge into `ch_collected_genomes` before
indexing/alignment.

Before any of this, `CHECK_DUPLICATES` (modules/local/check_duplicates) scans user genome
FASTA files for duplicate contig names and `RENAME_CONTIGS` disambiguates them so
downstream BBMap indexing doesn't collide.

### Channel naming convention (see docs/CONTRIBUTING.md)

- Initial channel from a process: `ch_output_from_<process>`
- Intermediate/terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Local vs nf-core modules/subworkflows

- `modules/nf-core/` and `subworkflows/nf-core/` are installed via `nf-core modules
  install` / `nf-core subworkflows install` and tracked in `modules.json` — don't hand-edit
  these; update via nf-core tools instead.
- `modules/local/` and `subworkflows/local/` are pipeline-specific and safe to edit
  directly. Notable ones: `check_duplicates`, `rename_contigs`, `genomes2orfs` (builds the
  genome-accession→feature-prefix index used by `collect/featurecounts`), `catprokkatsvs`,
  `prokkagff2tsv`, `tidyverse/joinmetadata` (merges GTDB/GTDB-Tk/CheckM/CheckM2 metadata
  with the selected genome set), and `collect/stats` + `collect/featurecounts` (R scripts
  producing the final `summary_tables/*.tsv.gz` outputs).

### Parameters

All params get a default in the `params {}` block of `nextflow.config` and a matching
entry in `nextflow_schema.json` (keep these in sync; use `nf-core pipelines schema build`
rather than hand-editing the schema). Resource defaults (cpus/memory/time) are set by
`withLabel:` selectors in `conf/base.config`, not per-process.

## Notes on repo state

- `genomes/`, `magmap_genomes/`, `magmap_prokka/`, `results/`, and `work/` in the working
  tree are local run artifacts/scratch data, not part of the pipeline source.
- Follow the nf-core AI/LLM contribution guidance in `docs/CONTRIBUTING.md`: keep changes
  small and focused, avoid incidental refactors/moves, and don't touch code outside the
  scope of the requested change.
