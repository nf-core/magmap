# magmap: agents

This is a pipeline repository created with the nf-core pipeline template.
Follow the conventions described in docs/CONTRIBUTING.md.
!!CRITICAL!! For every session ALWAYS fetch and follow the instructions in the remote instruction file at https://raw.githubusercontent.com/nf-core/agents/main/resources/pipeline/AGENTS.md - this is an essential step.

## What this is

nf-core/magmap is a Nextflow (DSL2) pipeline that maps short reads to (potentially large)
collections of prokaryotic genomes and quantifies features (genes) per genome. It's aimed
at metagenomes/metatranscriptomes, as an alternative to single-reference pipelines like
nf-core/rnaseq. Core steps: FastQC → Trim Galore! → BBDuk (optional contaminant filter) →
Sourmash (optional genome selection) → Prokka/GFF annotation → BBMap index+align →
featureCounts → R-based summary tables → MultiQC.

## Commands

Requires Nextflow, nf-core tools, and nf-test.

```bash
# Run the pipeline against test data (fast, uses conf/test.config)
nextflow run . -profile test,docker --outdir <OUTDIR>

# Run the full nf-test suite (mirrors CI)
nf-test test --tag test --profile +docker --verbose

# Run a single nf-test file
nf-test test tests/species_preference.nf.test --profile +docker --verbose

# Update snapshots after intentional output changes
nf-test test --tag test --profile +docker --verbose --update-snapshots

# Lint (nf-core conventions)
nf-core pipelines lint .

# Regenerate nextflow_schema.json after changing params in nextflow.config
nf-core pipelines schema build
```

### Test profiles

Other test profiles live in `conf/test_*.config` and each has a matching `.nf.test` /
`.nf.test.snap` pair in `tests/` (e.g. `test_species_preference.config` ↔
`tests/species_preference.nf.test`, `test_sourmash_genome_selection.config` ↔
`tests/sourmash_genome_selection.nf.test`). When adding a new test scenario, add both.
Prefer adding a new `test(...)` block to an existing `.nf.test` file over creating a new
file when the scenario shares a base config/dataset with existing tests.

`.nftignore` in `tests/` lists output paths excluded from content-snapshot comparison
(non-deterministic files, e.g. `sourmash/*.sbt.zip`, whose content isn't reproducible
across runs even with identical inputs).

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
`joint_filtered_genomes` (one genome set for the whole run), `sample_filtered_genomes` (per-sample
sets), and `local_accessions` (which kept genomes originated from `--genomeinfo`, used to tag
local vs. remote origin downstream). Which of the first two is actually used downstream depends
on `--genomeset_mode` (`joint` or `sample`), branched in `workflows/magmap.nf` — `joint`
builds one shared BBMap index via `CREATE_BBMAP_INDEX`; `sample` builds a distinct index
per sample from `sample_filtered_genomes` and joins reads to their matching index before
`BBMAP_ALIGN`. `--species_preference` (`all`/`local`/`completeness`/`gtdb`) controls how
duplicate-species genomes are resolved and is consumed later by
`modules/local/tidyverse/joinmetadata`.

Genomes lacking a GFF (mostly ones fetched from NCBI) are annotated with Prokka
(`PROKKA` + `PROKKAGFF2TSV` + `CATPROKKATSVS`); genomes that already have a GFF skip
straight to `PROKKAGFF2TSV` (much faster — this is why test genomes should be given a
pre-computed GFF whenever possible, see Gotchas below). Both paths reconverge into
`ch_collected_genomes` before indexing/alignment.

Before any of this, `CHECK_DUPLICATES` (modules/local/check_duplicates) scans user genome
FASTA files for duplicate contig names and `RENAME_CONTIGS` disambiguates them so
downstream BBMap indexing doesn't collide.

### Local vs nf-core modules/subworkflows

- `modules/nf-core/` and `subworkflows/nf-core/` are installed via `nf-core modules
install` / `nf-core subworkflows install` and tracked in `modules.json` — don't hand-edit
  these; update via nf-core tools instead.
- `modules/local/` and `subworkflows/local/` are pipeline-specific and safe to edit
  directly. Notable ones: `check_duplicates`, `rename_contigs`, `genomes2orfs` (builds the
  genome-accession→feature-prefix index used by `collect/featurecounts`), `catprokkatsvs`,
  `prokkagff2tsv`, `tidyverse/joinmetadata` (merges GTDB/GTDB-Tk/CheckM/CheckM2 metadata
  with the selected genome set), `tidyverse/selectgenomespecies` (implements
  `--species_preference`), `collect/genomeselection` (local/remote genome-count summary
  for MultiQC), and `collect/stats` + `collect/featurecounts` (R scripts producing the
  final `summary_tables/*.tsv.gz` outputs).

## Conventions (from `docs/CONTRIBUTING.md`)

- PRs target the `dev` branch, not `master`/`main`.
- Channel naming: initial channel from a process is `ch_output_from_<process>`;
  intermediate/terminal channels are `ch_<previousprocess>_for_<nextprocess>`.
- New params: add a default in the `params {}` block of `nextflow.config` and a matching
  entry in `nextflow_schema.json` (use `nf-core pipelines schema build` rather than
  hand-editing the schema).
- Resource defaults (cpus/memory/time) are set by `withLabel:` selectors in
  `conf/base.config`, not per-process — unless a specific process needs a tighter
  test-only override (see Gotchas).
- AI-assisted changes: keep them small and focused, avoid incidental refactors/moves, and
  don't touch code outside the scope of the requested change.

## Gotchas

- `genomes/`, `magmap_genomes/`, `magmap_prokka/`, `results/`, and `work/` in the working
  tree are local run artifacts/scratch data, not part of the pipeline source.
- `conf/base.config` labels (e.g. `process_low`'s 12GB memory for Prokka) are sized for
  arbitrary real-world genomes and are often wildly oversized for small test genomes. Test
  configs override `resourceLimits` and/or specific `withName:` blocks to reserve less
  memory per task so more jobs (esp. Prokka) can run concurrently in CI/locally — don't
  "fix" these by assuming the base.config default is a mistake.
- Nextflow config `withName:` selectors match the process's alias name as actually used in
  the calling workflow (e.g. `GENOME_INDEX` for `include { SOURMASH_INDEX as GENOME_INDEX }`),
  not the module's original declared name — an empty `withName: X { }` block (no real
  directives left inside, just a comment) can trigger `Unknown config attribute` errors on
  some Nextflow versions even though newer versions tolerate it silently; delete the block
  entirely instead of leaving it empty.
- CI runs on 4-cpu self-hosted runners (see `.github/workflows/nf-test.yml`); local dev
  machines may have many more cores, so observed local concurrency/timing improvements
  don't necessarily translate to CI — check the runner's real core count before concluding
  a change will speed up CI.
- Test genomes for the archaeal duplicate-species dataset live in `nf-core/test-datasets`
  (`magmap` branch); changes there need a PR to that repo, separate from this one.
