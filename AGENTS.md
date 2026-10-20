# magmap: agents

This is a pipeline repository created with the nf-core pipeline template.
Follow the conventions described in docs/CONTRIBUTING.md.
!!CRITICAL!! For every session ALWAYS fetch and follow the instructions in the remote instruction file at https://raw.githubusercontent.com/nf-core/agents/main/resources/pipeline/AGENTS.md - this is an essential step.

## Project notes

nf-core/magmap maps short reads to collections of prokaryotic genomes and quantifies features (genes) per genome.
It targets metagenomes and metatranscriptomes, as an alternative to single-reference pipelines such as nf-core/rnaseq.
The steps are FastQC, Trim Galore!, BBDuk (optional contaminant filter), Sourmash (optional genome selection), Prokka or Bakta annotation (for genomes without a GFF), BBMap index and align, featureCounts, summary tables and MultiQC.

### Commands

Needs Nextflow, nf-core tools and nf-test.

```bash
nextflow run . -profile test,docker --outdir <OUTDIR>
nf-test test --tag pipeline --profile=+docker --verbose
nf-test test tests/species_preference.nf.test --profile=+docker --verbose
nf-test test --tag pipeline --profile=+docker --verbose --update-snapshot
nf-core pipelines lint .
nf-core pipelines schema build
```

The commands are, in order: run the test profile, run the pipeline test suite as CI does, run one test file, update snapshots after an intended output change, lint, and rebuild `nextflow_schema.json` after changing parameters in `nextflow.config`.

### Entry points

`main.nf` runs `PIPELINE_INITIALISATION` (`subworkflows/local/utils_nfcore_magmap_pipeline`), the `NFCORE_MAGMAP` wrapper, `MAGMAP` (`workflows/magmap.nf`, the pipeline logic) and `PIPELINE_COMPLETION`.
Parameter validation, samplesheet parsing and the genome and index channels are set up in `PIPELINE_INITIALISATION`.
Look there first to see how a parameter becomes a channel.

### Genome selection

The set of reference genomes is built from three sources and then filtered per sample.

1. User genomes from `--genomeinfo` (a CSV with `accno,genome_fna,genome_gff`).
2. Local Sourmash indexes built from those genomes (skipped with `--skip_sourmash`).
3. Remote genomes found in NCBI assembly summary files (`--remote_genome_sources`) through user-supplied Sourmash indexes (`--indexes`), then fetched with `WGET`.

The logic is in `subworkflows/local/sourmash/main.nf` (workflow `SOURMASH`).
It emits `joint_filtered_genomes` (one set for the run), `sample_filtered_genomes` (one set per sample) and `local_accessions` (which kept genomes came from `--genomeinfo`).
`--genomeset_mode` (`joint` or `sample`) selects which set `workflows/magmap.nf` uses.
In `joint` mode `CREATE_BBMAP_INDEX` builds one shared index.
In `sample` mode it builds one index per sample and joins the reads to their index before `BBMAP_ALIGN`.
`--species_preference` (`all`, `local`, `completeness`, `gtdb`) controls how duplicate species are resolved, in `modules/local/tidyverse/selectgenomespecies`.
`--annotator` selects Prokka or Bakta.

Genomes without a GFF are annotated (`PROKKA` or `BAKTA_BAKTA`, then `PROKKAGFF2TSV` and `CATPROKKATSVS`).
Genomes with a GFF go straight to `PROKKAGFF2TSV`, which is much faster, so test genomes should come with a precomputed GFF.
Both paths join in `ch_collected_genomes` before indexing and alignment.
Before that, `CHECK_DUPLICATES` finds duplicate contig names in the user genomes and `RENAME_CONTIGS` makes them unique so that BBMap indexing does not collide.

### Local modules

- `duckdb/splitfeaturecounts` and `duckdb/joinfeaturecountsaccno`: split the featureCounts table per feature type and add the genome accession, in DuckDB so that memory stays bounded for large genome collections.
- `genomes2orfs`: builds the accession-to-ORF index.
- `tidyverse/joinmetadata`: merges GTDB, GTDB-Tk, CheckM and CheckM2 metadata with the selected genomes.
- `collect/genomeselection`: local and remote genome counts for MultiQC.
- `collect/unassignedcounts`: unassigned read counts for the summary tables.

`CUSTOM_COLLECTFEATURECOUNTS` and `CUSTOM_COLLECTSTATS` are nf-core modules that produce the final `summary_tables/*.tsv.gz`.

### Tests

Each `conf/test_*.config` has a matching `tests/*.nf.test` and `tests/*.nf.test.snap`.
Add a new `test(...)` block to an existing `.nf.test` file when the scenario shares its configuration and dataset.
`tests/.nftignore` lists output paths excluded from snapshot comparison because their content is not reproducible, for example `sourmash/*.sbt.zip`.
The pipeline snapshots list the files in `summary_tables` by name only, so a module test must check the table contents.

### Conventions

These come from `docs/CONTRIBUTING.md`.

- Pull requests target `dev`, not `master`.
- Name the first channel from a process `ch_output_from_<process>`, and the others `ch_<previousprocess>_for_<nextprocess>`.
- Give a new parameter a default in the `params {}` block of `nextflow.config` and an entry in `nextflow_schema.json`, made with `nf-core pipelines schema build` instead of by hand.
- Set resource defaults with `withLabel:` selectors in `conf/base.config`, not per process, unless a process needs a tighter test-only override.
- Keep AI-assisted changes small and focused, avoid incidental refactors and moves, and leave code outside the requested change alone.

### Gotchas

- `genomes/`, `magmap_genomes/`, `magmap_prokka/`, `results/` and `work/` in the working tree are local run artefacts, not pipeline source.
- The labels in `conf/base.config` are sized for real genomes and are far too large for test genomes. Test configs lower `resourceLimits` or set `withName:` blocks so that more tasks run at once. Do not change `base.config` to fix this.
- A `withName:` selector matches the alias used in the calling workflow, for example `GENOME_INDEX` for `include { SOURMASH_INDEX as GENOME_INDEX }`. An empty `withName:` block fails with `Unknown config attribute` on some Nextflow versions, so delete it.
- CI runs on 4-cpu runners. Local timings on a larger machine do not carry over.
- The test genomes live in the `magmap` branch of `nf-core/test-datasets`. Changes there need a separate pull request.
