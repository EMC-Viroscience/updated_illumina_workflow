# Illumina Workflow for Metagenomic Data Processing

**Adapted by:** Divyae Kishore Prasad  
**Original Workflow by:** Nathalie Worp and David Nieuwenhuisje  
**Development Period:** July–October 2024

## Overview

This repository is a Snakemake workflow for processing Illumina sequencing data in a metagenomic context. The pipeline includes steps for quality control, human read filtering, de novo assembly, annotation, and result summarization. The workflow is designed to be resource-aware, modular, and easy to configure, with outputs organized dynamically based on the current date.

## Features

- **Dynamic Output Folder:** Automatically creates an output folder named based on the current date (e.g., `processed_ddmmyy`).
- **Quality Control & Deduplication:** Uses [fastp](https://github.com/OpenGene/fastp) for QC and duplicate removal.
- **Human Read Filtering:** Removes contaminating human reads with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) and [samtools](http://www.htslib.org/).
- **De Novo Assembly:** Assembles filtered reads using [metaSPAdes](https://cab.spbu.ru/software/spades/).
- **Annotation:** Annotates assembled contigs with [DIAMOND BLASTX](https://github.com/bbuchfink/diamond) (using an e-value of `10-5`).
- **Mapping & Statistics:** Maps reads back to contigs and calculates coverage and mapping statistics using samtools and [seqkit](https://bioinf.shenwei.me/seqkit/).
- **Result Organization:** Merges results and creates renamed and linked annotation files for easier downstream analysis.
- **Rule Prioritization:** Uses rule priorities to help manage scheduling under limited resources.


## Installation

1. **Clone the Repository:**

    ```bash
    git clone https://github.com/divprasad/updated_illumina_workflow.git
    cd updated_illumina_workflow
    ```

2. **Set Up a Conda Environment (Recommended):**

    ```bash

    conda create -n illumina_workflow python=3.8 snakemake fastp bwa-mem2 samtools spades diamond seqkit -c bioconda -c conda-forge
    conda activate illumina_workflow

    ```

3. **Additional Dependencies:**

   The workflow uses optional custom scripts (e.g., `multiL_fasta_2singleL.py`).

## Configuration

The workflow uses a configuration variable `OUTPUT_FOLDER` (set in the Snakefile) to determine the output directory. If nothing is provided, then by default, it creates a folder named in the format `processed_ddmmyy` based on the current date. You can override this by providing a configuration file.

## Usage

**Check the `execute-and-log.sh` for easily launching the script with logging**

Alternatively, run the workflow by specifying the number of cores and the available memory. For example:

```bash

snakemake -s updated_illumina_workflow.smk \  # Specify the workflow file
    --config OUTPUT_FOLDER="processed_now" \  # Set the output directory
    --resources mem_gb=192 \  # Allocate 192GB of memory for the workflow
    # typically allocate mem_gb= 16x{cores}
    --cores 24 \  # Use up to 24 CPU cores
    --rerun-triggers mtime \  # Force rerun if input files are modified (based on modification time)
    --rerun-incomplete \  # Retry any incomplete jobs from previous runs
    --latency-wait 30 \  # Wait up to 30 seconds for input files (useful for network filesystems)
    --max-jobs-per-second 2 \  # Limit job submission rate to 2 per second
    --max-status-checks-per-second 4  # Limit status checks to 4 per second

```

##  Workflow Structure

1. **Raw Data Linking (softlink_raw):**
Soft-links raw FASTQ files from the raw_data/ directory into a dynamically generated output folder.

2. **Quality Control (QC_after_dedup):**
Runs fastp for QC and deduplication, outputting cleaned and deduplicated reads along with QC reports.

3. **Human Read Filtering (filter_human):**
Filters out human reads by aligning against the human genome reference.

4. **Assembly (assemble_filtered):**
Assembles filtered reads with metaSPAdes, followed by contig post-processing (e.g., reformatting FASTA files).

5. **Annotation (blastx_assembled):**
Annotates the assembled contigs using DIAMOND BLASTX.

6. **Mapping (map_reads_to_contigs):**
Maps reads back to the assembled contigs to generate BAM files.

7. **Statistics & Merging:** Generates coverage and read mapping statistics. Merges classification, coverage, and read mapping results into a final summary.

8. **File Organization (store_completed_annotation_files):** Renames and organizes the final annotation files and creates central links for downstream access.

## Customization

**Output Folder:**
You can provide a custom OUTPUT_FOLDER via a configuration file.

**Resource Allocation:**
The number of threads and memory for each rule are set dynamically. You can adjust these values using the threads and resources directives in each rule.

**Rule Priorities:**
Some rules have assigned priorities (e.g., assemble_filtered has priority 1, blastx_assembled has priority 2) to help guide the scheduler. You can adjust these if needed.

**Acknowledgements**
Nathalie Worp and David Nieuwenhuisje – for the original Illumina workflow that was further adapted here.
