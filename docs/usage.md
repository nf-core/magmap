# nf-core/magmap: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/magmap/usage](https://nf-co.re/magmap/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

Magmap is a workflow designed for mapping metatranscriptomic and metagenomic reads onto a group of genomes.
The collection of genomes can either be specified directly using a table (see the `--genomeinfo` parameter) or be the result of filtering with Sourmash.
The latter can use either the genomes specified by `--genomeinfo`, a "sketch index" pointing to genomes available for instance at NCBI (see the `--indexes` parameter) or a combination, to identify a smaller set to map to.
Genome files provided with `--genominfo` must include contigs in fasta format and optionally gff files (Prokka format).
Any genome for which a gff file is missing will be annotated with Prokka.
The pipeline can take output files from CheckM, CheckM2 and GTDB-Tk as input, and will provide processed output from these tools.
Note that the pipeline can map to any collection of genomes, including single genomes and isolates.

## Running the workflow

### Quickstart

A typical command for running the workflow is:

```bash
nextflow run nf-core/magmap -profile docker --outdir results/ --input samples.csv --genomeinfo localgenomes.csv
```

### Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It must be a comma-separated file with 3 columns, and a header row as shown in the examples below

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
T0a,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
T0b,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
T0c,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Genomes input

A second file input is the genome input sheet, which is specified with the option `--genomeinfo`. The file is a `.csv` file and it can contain three columns: `accno`, `genome_fna`, `genome_gff`. The first two are mandatory, while the third, `genome_gff`, is not.

```csv title="samplesheet.csv"
accno,genome_fna,genome_gff
GCA_002688505,./genomes/GCA_002688505.fna,./genomes/GCA_002688505.gff
GCA_002688515,./genomes/GCA_002688515.fna,
```

N.B. The pipeline assumes gff files have the same format as is output by Prokka.
Any genome used by the pipeline for which a gff file is not found will be annotated with Prokka to produce a gff file.

| Column       | Description                                                                                                                               |
| ------------ | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `accno`      | Accession number. For your local genomes, write the name of them.                                                                         |
| `genome_fna` | Full path to Fasta file that contains nucleotide sequences of your genome. File can be gzipped and have the extension ".fna.gz" or "fna". |
| `genome_gff` | Full path to gff file of your genome. File can be gzipped and have the extension ".gff.gz" or ".gff".                                     |

### Other inputs (optional)

Magmap can handle several types of input that can be used for different purposes.

#### Indexes input

In addition, or instead of, providing a genome file with genomes to map to, you can provide a [Sourmash](https://sourmash.readthedocs.io/en/latest/)
index file that points to genomes.
Sourmash will be run using the index files and matching genomes will be downloaded, annotated with Prokka and mapped to by the pipeline.
For this to work, entries in the Sourmash index need to point to NCBI assemblies with accessions in the format: `GC[A-Z]_[0-9]+\.[0-9]+`.
The indexes input is used by Sourmash to select genomes that can be downloaded in a second step and added to the pipeline.
It is provided with the `--indexes` parameter and it is a path (local or remote).

Particular examples of Sourmash index files are those prepared by the authors of Sourmash, which can be found
[here](https://sourmash.readthedocs.io/en/latest/databases.html).

```bash
nextflow run nf-core/magmap -profile docker --outdir results/ --input samples.csv --genomeinfo localgenomes.csv --indexes 'https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k21.sbt.zip'
```

N.B. More than one index file can be provided, separated by commas.

#### Genomes metadata input

Magmap accepts several types of metadata files that provides information about the genomes that you will use in the pipeline. Each type of metadata can be used to get different information about your genomes. At the moment, Magmap can handle output from CheckM/CheckM2 and gtdb-tk as well as standard GTDB metadata files. Magmap will merge all these tables and create a new one: each row will correspond to a genome (based on its accno) followed by several columns.

##### gtdb_metadata

With this parameter, you can supply a file like the GTDB metadata files provided on their official [website](https://gtdb.ecogenomic.org/), e.g. [`bac120_metadata_r220.tsv.gz`](https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz). You can either use their files directly or make a custom one. If you want make your own table, fill up the following columns: `accno`, `checkm_completeness`, `checkm_contamination`, `checkm_strain_heterogeneity`, `contig_count`, `genome_size`, `gtdb_genome_representative`,gtdb_representative`, `gtdb_taxonomy`.

##### gtdb-tk metadata

This file should be formatted like the GTDB-Tk output, [see](https://ecogenomics.github.io/GTDBTk/files/summary.tsv.html).

##### checkM/CheckM2 metadata

This file should be formatted like the checkM2 output, [see](https://github.com/nf-core/test-datasets/blob/magmap/testdata/checkm2.quality_report.tsv).

### Check duplicates

The pipeline will perform validation checks to see if there are any duplicate names among the genomes that the user provides. If there are duplicates, the pipeline will stop and return a file with the contig names that needs to be changed in their name in order to work. This is done to avoid overlapping in the following steps (e.g. same prokka output for the protein sequences and the gffs).

### Filter/remove sequences from the samples (e.g. rRNA sequences with SILVA database)

The pipeline can remove potential contaminants using the BBduk program.
Specify a fasta file, gzipped or not, with the `--sequence_filter <sequences>.fasta` parameter.
For further documentation, see the [BBduk official website](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

```bash
nextflow run nf-core/magmap -profile docker --outdir results/ --input samples.csv --genomeinfo localgenomes.csv --sequence_filter path/to/file
```

### Sourmash (optional)

With [Sourmash](https://sourmash.readthedocs.io/en/latest/index.html) you can filter the genomes to be used by magmap in the mapping step. This function is optional but can speed up the process and let you get a better genomes /reads mapping ratio since you are removing all the genomes that are not passing the threshold (that you can select).

```bash
nextflow run nf-core/magmap -profile docker --outdir results/ --input samples.csv --genomeinfo localgenomes.csv --sourmash true
```

### ORF caller option

The pipeline uses [Prokka](https://github.com/tseemann/prokka) to call genes/ORFs from the genomes. This is suitable for prokaryotes and it provides a gff as output for downstream analysis. It also performs functional annotation of ORFs.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/magmap --input ./samplesheet.csv --outdir ./results --genomeinfo ./genomes.csv -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/magmap -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/magmap
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/magmap releases page](https://github.com/nf-core/magmap/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organization are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
