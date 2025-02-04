# Illumina Workflow for Metagenomic Data Processing

**Adapted by:** Divyae Kishore Prasad  
**Original Workflow by:** Nathalie Worp and David Nieuwenhuisje  
**Development Period:** Jul'24–Feb'25

---

# Overview

This repository is a Snakemake workflow for processing Illumina sequencing data in a metagenomic context. The end-to-end workflow, `updated_illumina_workflow.smk` incorporates steps for going from raw reads to taxonomic annoation: quality control, human read filtering, de novo assembly, annotation, and result summarization. It is designed to be resource-aware, modular, and easy to configure, with outputs organized dynamically based on the current date.

---

# Workflow

1. **Dynamic Output Folder & Raw Data Linking**  
   - **Feature:** Automatically creates an output folder named based on the current date (e.g., `processed_ddmmyy`).  
   - **Step (softlink_raw):** Soft-links or copies raw FASTQ files from `raw_data/` into this dynamically generated folder for straightforward data management.

2. **Quality Control & Deduplication**  
   - **Feature:** Uses [fastp](https://github.com/OpenGene/fastp) to perform both quality trimming and deduplication in a single run.  
   - **Step (QC_after_dedup):** Outputs cleaned, deduplicated reads, and provides an HTML/JSON report detailing quality metrics.

3. **Human Read Filtering**  
   - **Feature:** Removes contaminating human reads by aligning against the human reference genome with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2), then uses [samtools](http://www.htslib.org/) to retain only unmapped reads.  
   - **Step (filter_human):** Produces filtered reads free of host contamination to improve downstream assembly and annotation accuracy.

4. **De Novo Assembly**  
   - **Feature:** Assembles filtered reads using [metaSPAdes](https://cab.spbu.ru/software/spades/), optimized for metagenomic data.  
   - **Step (assemble_filtered):** Generates assembled contigs; includes post-processing such as renaming contigs for clarity or converting multi-line FASTA to single-line FASTA.

5. **Annotation**  
   - **Feature:** Runs high-speed homology searches with [DIAMOND BLASTX](https://github.com/bbuchfink/diamond) (e-value of `10-5`) to annotate contigs against a protein database.  
   - **Step (blastx_assembled):** Produces a tab-delimited output with taxonomic and functional information for each contig.

6. **Mapping & Statistics**  
   - **Feature:** Maps reads back to the assembled contigs (with `bwa-mem2` + `samtools`) and uses [seqkit](https://bioinf.shenwei.me/seqkit/) for read-count statistics.  
   - **Steps (map_reads_to_contigs, Statistics & Merging):**  
     - Generates coverage information and BAM files,  
     - Merges coverage and annotation results into final summary tables for each sample.

7. **Result Organization**  
   - **Feature:** Creates renamed, centrally linked annotation and summary files for easier access and downstream analysis.  
   - **Step (store_completed_annotation_files):** Renames `completed_{sample}_annotation.tsv` files, links them in a central `annotations/` folder, and cleans up temporary files.

8. **Rule Prioritization**  
   - Certain rules (e.g., **blastx_assembled**, **assemble_filtered**) have assigned priorities to help the Snakemake scheduler manage tasks efficiently under limited resources.

---

# Installation and Quick Start

1. **Clone the Repository:**

    ```bash
    git clone https://github.com/divprasad/updated_illumina_workflow.git
    cd updated_illumina_workflow
    ```

2. **Set Up a Conda Environment (Recommended):**

    ```bash

    conda env create -f environment.yaml
    conda activate illumina_workflow

    ```

 3. **Set working directory:**
    Copy your raw Illumina paired end read data (FASTQ files) into the following structure:   the current folder in this format
    ```
    raw_data/{run}/{sample}_R1_001.fastq.gz
    raw_data/{run}/{sample}_R2_001.fastq.gz
    ```
    This is critical as wildcards are populated from raw_data {run} and {samples} that meet this format.

  4. **Minimal execution of the Workflow**  
     ```bash
     snakemake -s updated_illumina_workflow.smk \
         --cores 8
     ```
     This will start running the pipeline with 8 cores, producing results in a directory named `processed_ddmmyy` by default.

---

# Project Structure & Key Outputs

**Example project layout** after cloning the repo and executing the workflow:

```
updated_illumina_workflow/
├── raw_data/
│   └── runXYZ/
│       ├── sampleA_R1_001.fastq.gz
│       └── sampleA_R2_001.fastq.gz
├── updated_illumina_workflow.smk
├── environment.yaml
└── processed_ddmmyy/  # automatically created on running Snakemake
    ├── runXYZ/
    │   └── sampleA/
    │       ├── raw/
    │       ├── dedup_qc/
    │       ├── filtered/
    │       ├── assembly/
    │       ├── mappings/
    │       ├── summary/
    │       └── ...
    └── ...
```

**Key points**:
- **`raw_data/`** holds the initial FASTQ files.
- **`processed_ddmmyy/`** is generated by the workflow and contains processed outputs, organized by subfolders for each `{run}` and `{sample}`, such as `{OUTPUT_FOLDER}/{run}/{sample}`

## Output Directories & Files

1. **`dedup_qc/`**  
   - Contains the **QC**ed and **deduplicated** reads and fastp reports.

2. **`filtered/`**  
   - Holds the **human-filtered** reads after removing host contamination.

3. **`assembly/contigs.fasta`**  
   - Final **assembled** contigs from metaSPAdes.

4. **`{sample}_annotation.tsv`**  
   - DIAMOND BLASTX annotation results for each assembled contig.

5. **`mappings/{sample}_mappings.bam`**  
   - BAM files with reads mapped back to contigs.

6. **`summary/`**  
   - Summary tables of coverage, read statistics, and merged annotation data for quick reference.

---

# Advanced usage and configuration

**To run the workflow with logging, use** `execute-and-log.sh`, which automates the process.

Alternatively, you can launch the workflow by specifying the number of cores, memory, and other parameters. For example:

```bash
snakemake -s updated_illumina_workflow.smk \
    --resources mem_gb=192 \
    --cores 24 \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --latency-wait 30 \
    --max-jobs-per-second 2 \
    --max-status-checks-per-second 4
```
- **Default Output**: If no custom folder is specified, Snakemake will create an output directory named `processed_ddmmyy`.  
- **cores**: Use upto 24 CPU cores
- **resources mem_gb**: Allocate 192 GB memory (roughly 16× the number of cores).  
- **rerun-triggers**: Force a rerun if modification times indicate inputs have changed.  
- **rerun-incomplete**: Rerun incomplete jobs from previous executions.
- **latency-wait**: Wait up to 30 seconds for input files (useful for network filesystems)
- **Job submission and status check rate**: Limit new job submission rate (`max-jobs-per-second`) to 2 new jobs per second; limit job status check rate (`max-status-checks-per-second`) to 4 checks per second



## Config file
The workflow reads configuration variables from a user-provided config file or passing flags to Snakemake

- **Default**: If no configuration is passed, the workflow outputs to `processed_ddmmyy`.
- **Override**: You can override settings by passing flags to Snakemake or using a YAML config file.

For example, to specify a custom output folder and change the minimum contig length:

```bash
snakemake -s updated_illumina_workflow.smk \
    --config OUTPUT_FOLDER="processed_mydate" MIN_CONTIG_LEN="300" \
    --cores 16
```

> **Note**  
> If you wish to filter out shorter contigs, set `MIN_CONTIG_LEN` to your desired threshold. You may also need to uncomment lines related to `fil_renamed_contigs` in the **assemble_filtered** rule in the snakefile.

## Resource Allocation**
  - **Threads and Memory**: Each rule in the Snakefile can dynamically allocate threads and memory.  
  - **Adjusting Resources**: Modify the `threads:` or `resources:` directives inside each rule if you need more fine-grained control.  
  - **Cluster/Server Usage**: When running on HPC systems, you can specify job submission parameters (e.g., Slurm, PBS) using `--cluster` or a profile.

## Priorities
Some rules include **priority settings** to optimize execution order. This ensures that computationally intensive steps start earlier, preventing bottlenecks and reducing total runtime.

**Example Priority Assignments**  
- **`blastx_assembled`** → **Priority 2**  
  - As the most time-consuming step, it is executed first.  
- **`assemble_filtered`** → **Priority 1**  
  - This step runs soon after to avoid idle time.

> **Note**  
> By default, rules have priority 0. You can raise or lower priority levels as needed.

---

## **Acknowledgements**
Special thanks to **Nathalie Worp** and **David Nieuwenhuisje** for developing the original Illumina workflow, which was further adapted here.
