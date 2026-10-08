# nf-core/magmap roadmap

Last revised 2026-10-08.
Current release: 1.3.0. `dev` is 1.3.1dev and holds only #275 (AWS full-test Fusion setting).

## Ground rules

- Budget: fewer than 100 changed files per release (`git diff --stat master...dev | tail -1`).
- One theme per release. Small bug fixes may ride along when they touch files already in the diff.
- Eukaryote support is additive, so it ships as minor releases (1.4, 1.5, ...).
  nf-core semver makes only these breaking: renamed/removed params, changed _mandatory_ samplesheet columns, changed output formats or layout, removed features.
  A new optional `--genomeinfo` column is not breaking.
- 2.0 is reserved for the first batch of breaking changes.
- External gates (nf-core/modules, nf-core/test-datasets `magmap` branch) are listed as PR 0 and do not count against the budget.
- Tags: `#n` = GitHub issue or PR.

## Overview

| Release | Theme                                                                      | Type  | Est. files | Gated by                      |
| ------- | -------------------------------------------------------------------------- | ----- | ---------- | ----------------------------- |
| 1.4.0   | Eukaryotic transcriptomes from `--genomeinfo`, annotated with TransDecoder | minor | ~60-80     | test data, module install     |
| 1.5.0   | MMETSP as a remote source, with its quality metadata                       | minor | ~40-60     | #239 OSF link verified, 1.4.0 |
| 1.6.0   | Eukaryotic genomes: MetaEuk, EukCC, `species_preference`                   | minor | ~50-70     | 1.4.0, #241 open questions    |
| 2.0.0   | Breaking batch, starting with dropping Sourmash for user-provided genomes  | major | ~30-50     | -                             |

## 1.4.0 in depth

Goal: map reads against user-provided eukaryotic transcriptomes together with prokaryotic genomes, and count features per reference with the same summary tables.
Motivation: the metaT of a mixed 92-sample project mapped poorly because the eukaryotic fraction was invisible (#239).

### Contract

- `--genomeinfo` keeps `accno,genome_fna,genome_gff`. No mandatory column changes.
- A transcriptome without a GFF is annotated with TransDecoder (`transdecoder/longorf` + `transdecoder/predict`, already in nf-core/modules).
  Its `.transdecoder.gff3` has transcript-relative coordinates, which is the reference magmap maps to.
- A new optional `--genomeinfo` column, `sequence_type` (`genome` by default, or `transcriptome`), marks a row as a transcriptome.
- Everything downstream (BBMap, featureCounts, summary tables) is unchanged.

### Behaviour by input

| Row                | GFF given | Today                            | 1.4.0           |
| ------------------ | --------- | -------------------------------- | --------------- |
| prokaryotic genome | yes       | skip annotation                  | same            |
| prokaryotic genome | no        | Prokka or Bakta by `--annotator` | same            |
| transcriptome      | yes       | not supported                    | skip annotation |
| transcriptome      | no        | not supported                    | TransDecoder    |

### Invariants and traps found in the code

- `GENOMES2ORFS` extracts feature IDs with `grep -o 'ID=[A-Z_0-9]\+'`.
  TransDecoder IDs such as `ID=cds.<transcript>` contain lowercase letters and dots, so they would be cut or missed.
  Must be fixed and tested in 1.4.0, or counts will not join to genomes.
- The accession is taken from the file name with a `G.._[0-9.]+_` regex, and a name that does not match is used whole.
  Check that this equals the `--genomeinfo` `accno` for transcriptome files (#239 item 4 is probably smaller than the issue says).
- `--annotator` routing and `bakta_supported_only` use the GTDB domain.
  Transcriptomes must bypass it.
- Sourmash only touches user-provided genomes when `--skip_sourmash false` (default is true); remote selection via `--indexes` is separate.
  A transcriptome is never sketched or indexed and is always kept. With `--skip_sourmash false` it must be excluded from sketching (to verify in the code).
- `species_preference` and GTDB metadata are prokaryote-only. Transcriptomes bypass them in 1.4.0.
- A transcript is a short "genome" with one main feature, so confirm BBMap and featureCounts counts empirically on real test data (#239 item 5).

### Bug fixes and riders (touch the same files)

- #265: Apptainer in the container condition of local modules. `GENOMES2ORFS` is in the diff anyway.
- #263: report selected remote genomes missing from the NCBI summaries. Touches `subworkflows/local/sourmash`.
- #229: document annotation tool version drift, in `docs/usage.md`.
- #269: update the metro map with gffread, Bakta and TransDecoder. It belongs with the headline PR, because the map has to show the new step.

### PR plan

| PR  | Content                                                                                                                                                        | Est. files |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- |
| 0a  | nf-core/test-datasets `magmap` branch: 2-3 small eukaryotic transcriptomes, README entry                                                                       | external   |
| 0b  | nf-core/modules: check `transdecoder/*` are current                                                                                                            | external   |
| 1   | Riders: #265, #263, #229                                                                                                                                       | ~15        |
| 2   | `GENOMES2ORFS` ID regex generalisation, with test                                                                                                              | ~5         |
| 3   | Install transdecoder modules, new annotation branch in `workflows/magmap.nf`, `sequence_type` in `assets/schema_genomeinfo.json` and its parsing, params, docs | ~35        |
| 4   | Test profile `test_transcriptome` with `.nf.test` + snapshot, metro map (#269), CHANGELOG                                                                      | ~20        |

### Done when

- `test_transcriptome` profile and nf-test pass on both container engines.
- `nf-core pipelines lint` and `nextflow lint` pass at the minimum and latest Nextflow versions.
- `docs/usage.md` and `docs/output.md` describe transcriptome input and limitations (intronic and intergenic metaG reads do not map).
- The existing prokaryotic tests have unchanged snapshots.

## Later releases

### 1.5.0: MMETSP as a remote source

- New fetch mechanism: MMETSP is OSF/CAMERA-hosted, not in NCBI `assembly_summary.txt` (#239 item 1).
- Verify the OSF Sourmash index link before relying on it. The issue says it was found by search only.
- Quality metadata: reuse the published per-sample BUSCO matrix (`MMETSP_all_evaluation_matrix.csv`) instead of running anything (#241).
- Depends on 1.4.0.

### 1.6.0: eukaryotic genomes

- `metaeuk/easypredict` for genomes without annotation (#240). Check its `.gff` against `PROKKAGFF2TSV` first.
- EukCC as the CheckM analogue (#241). It is not in nf-core/modules, so that is a PR 0.
- `species_preference` for eukaryotic entries: reuse the modes or add a scheme (#241 open question).

## 2.0.0: breaking batch

- Drop Sourmash for user-provided genomes (sketching, local indexes, `--skip_sourmash`, `--sourmash_save_sourmash` for local genomes).
  Remote selection via `--indexes` stays.
  It is breaking per the nf-core spec (removal of a supported feature and params), so it cannot go into 1.x.
  `--skip_sourmash` already defaults to true, so most users are unaffected.
- Check which tests, docs and `species_preference` paths depend on local Sourmash before estimating (`sourmash_genome_selection`, `species_preference`).
- Collect other breaking items here as they come up, so users migrate once.

## Backlog

| #    | Item                                      | Note                                                       |
| ---- | ----------------------------------------- | ---------------------------------------------------------- |
| #184 | inStrain                                  | Modules exist, needs a design. Deferred until after 1.3.0. |
| #164 | eggnog-mapper and KOfamscan               | Separate theme (functional annotation, as in metatdenovo). |
| #87  | Skip steps for a concatenated genome file | Interacts with genome selection. Needs scoping.            |
| #34  | HMM profile scan                          | No description. Needs scoping.                             |

## Traceability

| Item                             | Release |
| -------------------------------- | ------- |
| #239 transcripts, local          | 1.4.0   |
| #239 MMETSP remote               | 1.5.0   |
| #240 TransDecoder                | 1.4.0   |
| #240 MetaEuk                     | 1.6.0   |
| #241 BUSCO metadata for MMETSP   | 1.5.0   |
| #241 EukCC, `species_preference` | 1.6.0   |
| #265, #263, #229, #269           | 1.4.0   |
| Drop local Sourmash              | 2.0.0   |
| #184, #164, #87, #34             | backlog |

## Decisions log

- 2026-10-08: eukaryote support ships as 1.4.0 (additive). 2.0 is kept for breaking changes. Milestone "Magmap 2.0" is to be renamed.
- 2026-10-08: transcriptomes come before genomes.
- 2026-10-08: budget is fewer than 100 files per release.
- 2026-10-08: transcriptomes are marked by an optional `--genomeinfo` column named `sequence_type` (`genome` default, `transcriptome`).
- 2026-10-08: Sourmash for user-provided genomes is dropped in 2.0, not 1.4 (removal is breaking). 1.4 only keeps transcriptomes out of it.
- 2026-10-08: milestone #6 renamed from "Magmap 2.0" to "Magmap 1.4".
- 2026-10-08: no separate agenda. The open issues are the list.
- 2026-10-08: small bug fixes ride along in 1.4.0 instead of a separate release.

## Open questions

- Does EukCC work on transcriptome assemblies, and what do the MMETSP BUSCO/Transrate scores map to in `species_preference`? (#241)
- Is the OSF Sourmash index for MMETSP real and loadable? (#239)
