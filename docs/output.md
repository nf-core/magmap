# nf-core/magmap: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and the results are organized as follow:

- [Summary tables](#summary-tables) - Tab separated tables ready for further analysis in tools like R and Python
- [Module output](#module-output)
  - [Preprocessing](#preprocessing)
    - [FastQC](#fastqc) - Read quality control
    - [Trim galore!](#trim-galore) - Primer trimming
    - [BBduk](#bbduk) - Filter out sequences from samples that matches sequences in a user-provided fasta file (optional)
  - [Filtering genomes](#filter-genomes) - Generate a list of genomes that will be used for the mapping
    - [Sourmash](#sourmash) - Output from Sourmash filtering of genomes.
  - [Prokka](#prokka) - Output from Prokka
  - [Genome fetching](#genome-fetching) - Genomes fetched from remote sources
  - [Quantification of genome features](#quantification-of-genome-features)
    - [BBmap](#bbmap) - Output from BBmap
    - [FeatureCounts](#featureCounts) - Output from FeatureCounts
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
- [MultiQC](#multiqc) - Aggregate report describing results

## Summary tables

Consistently named and formatted output tables in TSV format ready for further analysis.

<details markdown="1">
<summary>Output files</summary>

- `summary_tables/`
  - `magmap.overall_stats.tsv.gz`: Overall statistics from the pipeline, e.g. number of reads, number of called ORFs, number of reads mapping back to contigs/ORFs etc.
  - `magmap.<FEATURE>.counts.tsv.gz`: Read counts for `FEATURE` per ORF and sample.
  - `magmap.genome_metadata.tsv.gz`: Genome metadata from GTDB, GTDB-Tk and CheckM/CheckM2 if provided by the user.
  - `magmap.genomes2orfs.tsv.gz`: Translation table from ORF identifiers to genome identifiers.
  - `magmap.prokka-annotations.tsv.gz`: Annotation details extracted from GFF files.

</details>

## Module output

### Preprocessing

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/). FastQC is run as part of Trim galore! therefore its output can be found in Trim galore's folder.

<details markdown="1">
<summary>Output files</summary>

- `trimgalore/fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.

</details>

#### Trim galore!

[Trim galore!](https://github.com/FelixKrueger/TrimGalore) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primers.

<details markdown="1">
<summary>Output files</summary>

- `trimgalore/`: directory containing log files with retained reads, trimming percentage, etc. for each sample.
  - `*trimming_report.txt`: report of read numbers that pass trimgalore.

</details>

#### BBduk

[BBduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/) is a filtering tool that removes specific sequences from the samples using a reference fasta file.
BBduk is built-in tool from BBmap.

<details markdown="1">
<summary>Output files</summary>

- `bbmap/`
  - `*.bbduk.log`: a text file with the results from BBduk analysis. Number of filtered reads can be seen in this log.

</details>

### Filtering genomes

The Sourmash program can be used to prefilter genomes so that only genomes likely to be represented among the reads are passed to mapping.
In addition, Sourmash can be used to fetch remote genomes, see [usage docs](https://nf-co.re/magmap/usage#genome-input).
No output from Sourmash is enabled by default; the output is only used to select genomes for further processing.
Use [`--sourmash_save_sourmash`](https://nf-co.re/magmap/parameters/#sourmash_save_sourmash) to copy output files.

<details markdown="1">
<summary>Output files</summary>

- `sourmash/`
  - `*`: Output from Sourmash

</details>

### Prokka

[Prokka](https://github.com/tseemann/prokka) will be used to identify ORFs in any genomes for which a gff file is not provided.
In addition to calling ORFs (done with Prodigal) Prokka will functionally annotate the ORFs.
To make it easier to reuse already annotated genomes in other projects, output from Prokka is directed to subdirectories of the directory specified with the [`--prokka_store_dir` parameter](https://nf-co.re/magmap/parameters/#prokka_store_dir) (by default `prokka` in the working directory for the pipeline run).
Genomes already found in the directory specified, will be skipped by the Prokka step.

<details markdown="1">
<summary>Output files</summary>

- `prokka/`
  - `<accno>`
    - `*.ffn`: nucleotide fasta file output
    - `*.faa`: amino acid fasta file output
    - `*.gff`: genome feature file output

</details>

### Genome fetching

When the pipeline is run with [`--skip_sourmash false`](https://nf-co.re/magmap/parameters/#skip_sourmash) and one or more index files passed to [`--indexes`](https://nf-co.re/magmap/parameters/#indexes), remote genomes will be identified and downloaded to the directory specified by [`--genome_store_dir`](https://nf-co.re/magmap/parameters/#genome_store_dir) (by default `genomes` in the working directory for the pipeline run).

### Quantification of genome features

#### BBmap

Only logs are saved by default from the BBmap step.
To save the `.bam` files, use `--bbmap_save_bam` and to save the index, use `--bbmap_save_index`.

<details markdown="1">
<summary>Output files</summary>

- `bbmap/`
  - `bam/`
    - `<SAMPLE>.bam`: bam file for `SAMPLE`
  - `logs/`:
    - `<SAMPLE>.bbmap.log`: BBmap log for `SAMPLE`

<details markdown="1">
<summary>Output files</summary>

#### FeatureCounts

<details markdown="1">
<summary>Output files</summary>

- `featurecounts/`
  - `<SAMPLE>.<FEATURE>.featureCounts.tsv`: Counts for `SAMPLE` and `FEATURE`
  - `<SAMPLE>.<FEATURE>.featureCounts.tsv.summary`: Summary of counts for `SAMPLE` and `FEATURE`

<details markdown="1">

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::
