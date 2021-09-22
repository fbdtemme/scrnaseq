# nf-core/scrnaseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [nf-core/scrnaseq: Output](#nf-corescrnaseq-output)
    * [:warning: Please read this documentation on the [nf-core website](https://nf-co.re/scrnaseq/output)] (#warning-please-read-this-documentation-on-the-nf-core-website-httpsnf-corescrnaseqoutput)
    * [Introduction](#introduction)
    * [Pipeline overview](#pipeline-overview)
    * [FastQC](#fastqc-results)
    * [Kallisto & Bustools](#kallisto--bustools-results)
    * [STARsolo](#starsolo)
    * [Salmon Alevin & AlevinQC](#salmon-alevin--alevinqc)
    * [Salmon Alevin-fry](#salmon-alevin-fry)
    * [CellRanger](#cellranger)
    * [Other output data](#other-output-data)
    * [MultiQC](#multiqc)
    * [Pipeline information](#pipeline-information)

## FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain adapter sequence and regions of low quality.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

## Kallisto & Bustools

See [Kallisto](https://pachterlab.github.io/kallisto/about) for details about Kallisto and [Bustools](https://bustools.github.io/) for more information on BusTools.

The pipeline can analyze data from single cell rnaseq experiments and generates a set of folders with respective outputs from various steps of the analysis. For a detailed summary what the pipeline does specifically, please follow the [excellent tutorial](https://www.kallistobus.tools/getting_started.html) that also describes specific steps for downstream analysis of the generated matrices.

**Output directory: `results/kallisto`**

* `sample.count`
    * Contains the output from Kallisto bustools

* `sample_matrices`
    * Contains the matrix of counts in different formants and some QC plots

**Output directory: `results/index`**

* `kallisto`
    * Contains the index of the supplied (genome/transcriptome) fasta file

## STARsolo

**Output directory: `results/star`**

* `sample.Solo.out`
    * Contains the output from STARSolo

* `sample_matrices`
    * Contains the matrix of counts in different formants and some QC plots

**Output directory: `results/index`**

* `star`
    * Contains the index of the supplied genome fasta file

## Salmon Alevin & AlevinQC

**Output directory: `results/salmon`**

* `sample_alevin_results`
    * Contains the output from Salmon Alevin
* `sample_alevin_matrices`
    * Contains the matrix of counts in different formants and some QC plots
* `alevinqc`
    * Contains the QC web report generated from the Salmon Alevin output

**Output directory: `results/index`**

* `salmon`
    * Contains the index reference for Salmon Alevin

## Salmon Alevin-fry

**Output directory: `results/alevinfry`**

* `sample_quant_results`
    * Contains the output from Alevin-fry
* `sample_quant_matrices`
    * Contains the matrix of counts in different formants and some QC plots

**Output directory: `results/index`**

* `alevinfry`
    * Contains the index reference for Alevin-fry

## CellRanger

**Output directory: `results/cellranger`**

* `sample/outs`
    * Contains output from CellRanger
* `matrices`
    * Contains the matrix of counts in different formants and some QC plots

**Output directory: `results/index`**

* `cellranger`
    * Contains the index reference for CellRanger

## Other output data
  
## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
    * Documentation for interpretation of results in HTML format: `results_description.html`.
