# nf-core/magmap: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev - [YYYY-mm-dd]

### `Added`

- [#217](https://github.com/nf-core/magmap/pull/217) - Add `--annotator` parameter to allow Bakta as an alternative to Prokka for genomes lacking a gff (`prokka`, `bakta_supported_only` or `bakta_all`), with `--bakta_db` and `--bakta_store_dir` controlling the Bakta database and annotation output locations (by @erikrikarddaniel).
- [#212](https://github.com/nf-core/magmap/pull/212) - Add `--species_preference` parameter to allow choice of genome representing a species (by @erikrikarddaniel).
- [#212](https://github.com/nf-core/magmap/pull/212) - Do not publish Sourmash signature files (`*.sig` and `*.sig.zip`) even with `--sourmash_save_sourmash` (by @erikrikarddaniel).
- [#212](https://github.com/nf-core/magmap/pull/212) - Add summary of selected genomes to MultiQC report (by @erikrikarddaniel).
- [#219](https://github.com/nf-core/magmap/pull/219) - Add `--save_parquet` to also write the summary tables in Parquet format, alongside the default gzipped TSV (by @erikrikarddaniel).

### `Changed`

- [#222](https://github.com/nf-core/magmap/pull/222) - Define TPM calculation in docs/output.md (by @erikrikarddaniel).
- [#223](https://github.com/nf-core/magmap/pull/223) - Template update for nf-core/tools version 4.1.0, and update nf-core modules and subworkflows (by @erikrikarddaniel).
- [#215](https://github.com/nf-core/magmap/pull/215) - Template update for nf-core/tools version 4.0.3, and update nf-core modules and subworkflows (by @erikrikarddaniel).

### `Fixed`

- [#222](https://github.com/nf-core/magmap/pull/222) - Fix brittle R code in `COLLECT_STATS`/`COLLECT_GENOMESELECTION` ([#218](https://github.com/nf-core/magmap/issues/218), by @erikrikarddaniel).
- [#222](https://github.com/nf-core/magmap/pull/222) - Fix an intermittent `ConcurrentModificationException` in `COLLECT_GENOMESELECTION` caused by two channel subscribers sharing a mutable genome-accession list (by @erikrikarddaniel).
- [#217](https://github.com/nf-core/magmap/pull/217) - Report the Prokka version in the collated `versions.yml`/MultiQC report; it was silently never being collected (by @erikrikarddaniel).
- [#215](https://github.com/nf-core/magmap/pull/215) - Fix a CI issue affecting Singularity/Apptainer test jobs (by @erikrikarddaniel).
- [#223](https://github.com/nf-core/magmap/pull/223) - Keep the locally-reported Prokka version (`PROKKA_VERSION`) in sync with the vendored module's actual version, and reuse its container instead of a near-duplicate one (by @erikrikarddaniel).

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| MultiQC    | 1.34             | 1.35        |
| Samtools   | 1.23.1           | 1.24        |
| trimgalore | 2.1.0            | 2.3.0       |
| Nextflow   | 25.10.4          | 26.04.0     |
| Prokka     | 1.14.6           | 1.15.6      |

### `Deprecated`

## v1.1.0 - [2026-05-28]

### `Added`

- [#199](https://github.com/nf-core/magmap/pull/199) - update template v4.0.2, created with the [nf-core](https://nf-co.re/) template (by @danilodileo).
- [#189](https://github.com/nf-core/magmap/pull/189) - Add genome `accno` field to counts tables (by @erikrikarddaniel).
- [#188](https://github.com/nf-core/magmap/pull/188) - Add the possibility to run mapping with sample-specific genome sets (by @erikrikarddaniel).

### `Fixed`

- [#210](https://github.com/nf-core/magmap/pull/210) - Update nf-schema and remove genomes param warning. (by @erikrikarddaniel)
- [#201](https://github.com/nf-core/magmap/pull/201) - Changed name inside featureCounts output in order to obtain multiQC output for each feature (issue [#200](https://github.com/nf-core/magmap/issue/200)
- [#198](https://github.com/nf-core/magmap/pull/198) - Update the --indexes parameter so it takes a list of indexes (issue [#200](https://github.com/nf-core/magmap/pull/201); by @danilodileo).
- [#197](https://github.com/nf-core/magmap/pull/197) - Fix bug in sample-specific genome set mode (by @erikrikarddaniel).
- [#186](https://github.com/nf-core/magmap/pull/186) - Update pipeline to use topic channels for tool versions (by @erikrikarddaniel).
- [#186](https://github.com/nf-core/magmap/pull/186) - Make sure `nextflow lint` works for the pipeline (by @erikrikarddaniel).
- [#186](https://github.com/nf-core/magmap/pull/186) - Update pipeline to nf-core template v.3.5.2 (by @erikrikarddaniel).
- [#178](https://github.com/nf-core/magmap/pull/178) - Make sure output from samtools is saved in `<outdir>`.

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| MultiQC    | 1.32             | 1.33        |
| Subread    | 2.0.6            | 2.1.1       |
| trimgalore | 0.6.10           | 2.1.0       |
| Samtools   | 1.22.1           | 1.23.1      |

### `Deprecated`

## v1.0.0 - [2025-11-28]

Initial release of nf-core/magmap, created with the [nf-core](https://nf-co.re/) template.
