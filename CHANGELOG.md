# nf-core/magmap: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [2026-05-28]

### `Added`

- [#199](https://github.com/nf-core/magmap/pull/199) - update template v4.0.2, created with the [nf-core](https://nf-co.re/) template. (by @danilodileo)
- [#189](https://github.com/nf-core/magmap/pull/189) - Add genome `accno` field to counts tables (by @erikrikarddaniel).
- [#188](https://github.com/nf-core/magmap/pull/188) - Add the possibility to run mapping with sample-specific genome sets (by @erikrikarddaniel).

### `Fixed`

- [#201](https://github.com/nf-core/magmap/pull/201) - Changed name inside featureCounts output in order to obtain multiQC output for each feature (issue [#200](https://github.com/nf-core/magmap/issue/200)
- [#198](https://github.com/nf-core/magmap/pull/198) - Update the --indexes parameter so it takes a list of indexes (issue [#200](https://github.com/nf-core/magmap/pull/201); by @danilodileo).
- [#197](https://github.com/nf-core/magmap/pull/197) - Fix bug in sample-specific genome set mode (by @erikrikarddaniel).
- [#186](https://github.com/nf-core/magmap/pull/186) - Update pipeline to use topic channels for tool versions (by @erikrikarddaniel).
- [#186](https://github.com/nf-core/magmap/pull/186) - Make sure `nextflow lint` works for the pipeline (by @erikrikarddaniel).
- [#186](https://github.com/nf-core/magmap/pull/186) - Update pipeline to nf-core template v.3.5.2 (by @erikrikarddaniel).
- [#178](https://github.com/nf-core/magmap/pull/178) - Make sure output from samtools is saved in `<outdir>`.

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- | ------ |
| MultiQC | 1.32             | 1.33        |
| +       | Subread          | 2.0.6       | 2.1.1  |
| +       | trimgalore       | 0.6.10      | 2.1.0  |
| +       | Samtools         | 1.22.1      | 1.23.1 |

### `Deprecated`

## v1.0.0 - [2025-11-28]

Initial release of nf-core/magmap, created with the [nf-core](https://nf-co.re/) template.
