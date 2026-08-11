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
  - [Bakta](#bakta) - Output from Bakta
  - [Genome fetching](#genome-fetching) - Genomes fetched from remote sources
  - [Quantification of genome features](#quantification-of-genome-features)
    - [BBmap](#bbmap) - Output from BBmap
    - [FeatureCounts](#featurecounts) - Output from FeatureCounts
    - [Samtools](#samtools) - Output from Samtools
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
- [MultiQC](#multiqc) - Aggregate report describing results

## Summary tables

Consistently named and formatted output tables in TSV format ready for further analysis.
With [`--save_parquet`](https://nf-co.re/magmap/parameters/#save_parquet), the same tables are also written as [Parquet](https://parquet.apache.org/) files, alongside the TSVs.

<details markdown="1">
<summary>Output files</summary>

- `summary_tables/`
  - `magmap.overall_stats.tsv.gz`: Overall statistics from the pipeline, e.g. number of reads, number of called ORFs, number of reads mapping back to contigs/ORFs etc.
  - `magmap.<FEATURE>.counts.tsv.gz`: Read counts and TPMs for `FEATURE` per ORF and sample. TPMs are calculated as: `r = count/length; tpm = r/sum(r)` over each sample, i.e. a length-corrected relative abundance.
  - `magmap.genome_metadata.tsv.gz`: Genome metadata from GTDB, GTDB-Tk and CheckM/CheckM2 if provided by the user.
  - `magmap.genome_selection.tsv.gz`: Per genome set (one per sample when `--genomeset_mode sample` is used, a single set otherwise), which genomes were selected and whether each originated from `--genomeinfo` (local) or was fetched from NCBI (remote).
  - `magmap.genomes2orfs.tsv.gz`: Translation table from ORF identifiers to genome identifiers.
  - `magmap.prokka-annotations.tsv.gz`: Annotation details extracted from GFF files.
  - `*.parquet`: with `--save_parquet`, a Parquet copy of each of the above (same basename, `.parquet` extension instead of `.tsv.gz`).

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
No output from Sourmash is saved to `<outdir>` by default; the output is only used to select genomes for further processing.
Use [`--sourmash_save_sourmash`](https://nf-co.re/magmap/parameters/#sourmash_save_sourmash) to copy the `*.csv.gz` and `*.sbt.zip` output files (`*.sig` and `*.sig.zip` are not saved under `<outdir>` even with this parameter).

<details markdown="1">
<summary>Output files</summary>

- `sourmash/`
  - `*.csv.gz`: Comma-separated file with Sourmash statistics.

</details>

### Prokka

[Prokka](https://github.com/tseemann/prokka) identifies ORFs in genomes for which a gff file is not provided and [`--annotator`](https://nf-co.re/magmap/parameters/#annotator) routes to it -- see [Bakta](#bakta) below for the alternative.
In addition to calling ORFs (done with Prodigal) Prokka will functionally annotate the ORFs.
To make it easier to reuse already annotated genomes in other projects, output from Prokka is directed to subdirectories of the directory specified with the [`--prokka_store_dir` parameter](https://nf-co.re/magmap/parameters/#prokka_store_dir) (by default `magmap_prokka` in the working directory for the pipeline run).
Genomes already found in the directory specified, will be skipped by the Prokka step.

<details markdown="1">
<summary>Output files</summary>

- `magmap_prokka/`
  - `<accno>/`
    - `*.ffn.gz`: nucleotide fasta file output
    - `*.faa.gz`: amino acid fasta file output
    - `*.gff.gz`: genome feature file output

</details>

### Bakta

[Bakta](https://github.com/oschwengers/bakta) is an alternative to Prokka for genomes lacking a gff, selected with [`--annotator bakta_supported_only` or `--annotator bakta_all`](https://nf-co.re/magmap/parameters/#annotator) -- see [Choosing an annotator: Prokka or Bakta](https://nf-co.re/magmap/usage#choosing-an-annotator-prokka-or-bakta) in the usage docs.
Bakta is only designed to annotate Bacteria; `bakta_supported_only` keeps Archaea (and genomes it can't classify by domain) on Prokka, while `bakta_all` sends everything to Bakta regardless.
As with Prokka, output is directed to subdirectories of the directory specified with [`--bakta_store_dir`](https://nf-co.re/magmap/parameters/#bakta_store_dir) (by default `magmap_bakta`), and genomes already found there are skipped.
The Bakta database itself is downloaded once into [`--bakta_db`](https://nf-co.re/magmap/parameters/#bakta_db) (by default `magmap_bakta_db`) and reused on subsequent runs.

<details markdown="1">
<summary>Output files</summary>

- `magmap_bakta/`
  - `<accno>/`
    - `*.ffn.gz`: nucleotide fasta file output
    - `*.faa.gz`: amino acid fasta file output
    - `*.gff3.gz`: genome feature file output
    - `*.tsv.gz`, `*.txt.gz`, `*.json.gz`, `*.embl.gz`, `*.gbff.gz`, `*.hypotheticals.*.gz`: Bakta's other native output formats

</details>

:::warning
`--prokka_store_dir`/`--bakta_store_dir` are reused across pipeline runs -- genomes already found there are skipped rather than re-annotated (see above).
This means genomes accumulated in a store directory over time may have been annotated with _different_ versions of Prokka or Bakta, even though only a single, current version is reported in the collated `versions.yml`/MultiQC report.
There is currently no automated check or report of per-genome annotation tool versions across a genome library (see [issue #229](https://github.com/nf-core/magmap/issues/229)).
If consistent annotation-tool versions across your whole genome collection matters, either start from an empty store directory, or track outside the pipeline which version annotated which genome.
:::

### Genome fetching

When the pipeline is run with [`--skip_sourmash false`](https://nf-co.re/magmap/parameters/#skip_sourmash) and one or more index files passed to [`--indexes`](https://nf-co.re/magmap/parameters/#indexes), remote genomes will be identified and downloaded to the directory specified by [`--genome_store_dir`](https://nf-co.re/magmap/parameters/#genome_store_dir) (by default `genomes` in the working directory for the pipeline run).
The exact behaviour is affected also by [--genomeset_mode](https://nf-co.re/magmap/parameters/#genomeset_mode) and [--species_preference](https://nf-co.re/magmap/parameters/#species_preference).

<details markdown="1">
<summary>Output files</summary>

- `bbmap/`
  - `*.genomes.txt`: Selected genomes, either for all samples or per sample depending on the `--genomeset_mode` parameter setting.

</details>

### Quantification of genome features

#### BBmap

[BBMap](https://sourceforge.net/projects/bbmap/) is a splice-aware global aligner for DNA and RNA sequencing reads. It provides detailed alignment statistics including mapped read counts, insert size distribution, and error rates. For further reading and documentation see the [BBMap documentation](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/).

Only logs and a manifest of the genome accessions included in each index are saved by default from the BBmap step.
To save the `.bam` files, use `--bbmap_save_bam`; to save the index, use `--bbmap_save_index`.

<details markdown="1">
<summary>Output files</summary>

- `bbmap/`
  - `bam/`
    - `<SAMPLE>.bam`: bam file for `SAMPLE`
  - `logs/`
    - `<SAMPLE>.bbmap.log`: BBmap log for `SAMPLE`
  - `<SAMPLE or "all">.genomes.txt`: genome accessions included in the BBMap index used for `SAMPLE` (`--genomeset_mode sample`) or for the whole run (`all`, `--genomeset_mode joint`)

</details>

#### FeatureCounts

[featureCounts](https://subread.sourceforge.net/) is a highly efficient read summarization program that counts mapped reads for genomic features such as genes and transcripts. It provides read count matrices that serve as input for downstream differential expression analysis. For further reading and documentation see the [featureCounts documentation](https://subread.sourceforge.net/featureCounts.html).

<details markdown="1">
<summary>Output files</summary>

- `featurecounts/`
  - `<SAMPLE>.<FEATURE>.featureCounts.tsv`: Counts for `SAMPLE` and `FEATURE`
  - `<SAMPLE>.<FEATURE>.featureCounts.tsv.summary`: Summary of counts for `SAMPLE` and `FEATURE`

</details>

#### Samtools

[Samtools](https://www.htslib.org) is a suite of programs for interacting with high-throughput sequencing data. It provides statistics on alignment quality, coverage, and read duplication. For further reading and documentation see the [Samtools documentation](http://www.htslib.org/doc/). Samtools statistics are provided alongside the alignment files in the BBMap output folder.

<details markdown="1">
<summary>Output files</summary>

- `samtools/`
  - `<SAMPLE>.flagstat`: Flagstat statistics for `SAMPLE`
  - `<SAMPLE>.idxstats`: Idxstats statistics for `SAMPLE`
  - `<SAMPLE>.stats`: Stats statistics for `SAMPLE`

</details>

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

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

The report also includes a "Genome selection" table with the number of local (user-provided) and remote (NCBI-fetched) genomes selected per genome set -- one row per sample when `--genomeset_mode sample` is used, a single row for the whole run otherwise. See `summary_tables/magmap.genome_selection.tsv.gz` for the full per-genome breakdown.
