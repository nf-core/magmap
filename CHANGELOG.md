# nf-core/magmap: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev - [YYYY-mm-dd]

### `Added`

- [#212](https://github.com/nf-core/magmap/pull/212) - Add `--species_preference` parameter to allow choice of genome representing a species (by @erikrikarddaniel).
- [#212](https://github.com/nf-core/magmap/pull/212) - Do not publish Sourmash signature files (`*.sig` and `*.sig.zip`) even with `--sourmash_save_sourmash` (by @erikrikarddaniel).
- [#212](https://github.com/nf-core/magmap/pull/212) - Add summary of selected genomes to MultiQC report (by @erikrikarddaniel).
- [#NN](https://github.com/nf-core/magmap/pull/NN) - Add `--annotator` parameter to allow Bakta as an alternative to Prokka for genomes lacking a gff (`prokka`, `bakta_supported_only` or `bakta_all`), with `--bakta_db` and `--bakta_store_dir` controlling the Bakta database and annotation output locations (by @erikrikarddaniel).

### `Changed`

- [#215](https://github.com/nf-core/magmap/pull/215) - Template update for nf-core/tools version 4.0.3, and update nf-core modules and subworkflows (by @erikrikarddaniel).

### `Fixed`

- [#NN](https://github.com/nf-core/magmap/pull/NN) - Report the Prokka version in the collated `versions.yml`/MultiQC report; it was silently never being collected (by @erikrikarddaniel).
- [#215](https://github.com/nf-core/magmap/pull/215) - Revert the `eWaterCycle/setup-apptainer` digest bumped in by the nf-core/tools 4.0.3 template sync -- that update is faulty and causes a namespace error in some pipelines' Singularity/Apptainer CI jobs (by @erikrikarddaniel).

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| MultiQC    | 1.34             | 1.35        |
| Samtools   | 1.23.1           | 1.24        |
| trimgalore | 2.1.0            | 2.3.0       |
| Nextflow   | 25.10.4          | 26.04.0     |

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
