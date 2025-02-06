# Illumina Workflow for Metagenomic Data Processing  


This repository contains a `snakemake` workflow for processing Illumina sequencing data, optimized and validated for metagenomics. The end-to-end workflow, `up_illumina_wf.snakefile`, includes steps to process raw reads into taxonomic annotation: quality control, human read filtering, de novo assembly, annotation, and result summarization. The Snakefile is designed to be resource-aware, modular, and easy to configure, with outputs dynamically organized based on the current date.  


## Workflow breakdown: from raw reads to annotation

1. **Dynamic output folder & raw data linking**  
   - Automatically creates an output folder named based on the current date (e.g., `processed_ddmmyy`).  
   - `hardlink_raw`: Hard-links or copies raw FASTQ files from `raw_data/`. The dynamically generated folder ensures easy data management.

2. **Quality control & deduplication**  
   - Performs adapter detection, quality trimming, and read deduplication in a single step using [fastp](https://github.com/OpenGene/fastp).  
   - `QC_after_dedup`: Generates cleaned, deduplicated reads, and provides an HTML/JSON report detailing QC metrics.

3. **Human read filtering**  
   - Removes (background contaminating) human reads by aligning against the human reference genome with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2), then uses [samtools](http://www.htslib.org/) to retain only unmapped reads.  
   - `filter_human`: Generates filtered reads free of host contamination to improve downstream assembly and annotation accuracy.

4. **De novo assembly**  
   - Assembles filtered reads using [metaSPAdes](https://cab.spbu.ru/software/spades/), optimized for metagenomic data.  
   - `assemble_filtered`: Generates assembled contigs; includes post-processing such as renaming contigs to include run and barcode information for clarity.

5. **Annotation**  
   - Runs high-speed homology searches with [DIAMOND BLASTX](https://github.com/bbuchfink/diamond) (e-value of `10-5`) to annotate contigs against a protein database.  
   - `blastx_assembled`: Produces a tab-delimited output with taxonomic and functional information for each contig.

6. **Mapping & statistics**  
   - Maps reads back to the assembled contigs (with `bwa-mem2` + `samtools`) and uses [seqkit](https://bioinf.shenwei.me/seqkit/) for read-count statistics.  
   - `map_reads_to_contigs`: Generates BAM files, creates coverage files, and mapping statistics.

7. **Merge results and organization**  
   - `merge_results`: Merges coverage and annotation results into summary tables.
   - `store_completed_annotation_files`: Creates renamed, centrally linked annotation for easier downstream analysis. It renames `completed_{sample}_annotation.tsv` files, links them in a central `annotations/` folder, and removes temporary files.

8. **Rule prioritization**  
   - Certain rules, such as `blastx_assembled` and `assemble_filtered`, have assigned priorities to optimize scheduling and execution.  


## Installation and quick start

1. **Clone the repository:**

    ```
    git clone https://github.com/divprasad/updated_illumina_workflow.git
    cd updated_illumina_workflow
    ```

2. **Set up a `conda` environment:**

    ```
    conda env create -f environment.yaml
    conda activate illumina_wf_env
    ```

3. **Set working directory:**
    Place raw Illumina paired-end FASTQ files in the following directory structure:
    ```
    raw_data/{run}/{sample}_R1_001.fastq.gz
    raw_data/{run}/{sample}_R2_001.fastq.gz
    ```
    This structure is critical for the workflow to recognize `run` and `sample` wildcards correctly. Ensure the file names follow this format.

4. **Run the workflow wrapper:**

    Launch `execute-and-log.sh`, which provides step-by-step guidance on executing `snakemake` and manages execution logging.  

5. Alternatively, **launch the workflow with minimal settings:**  
    ```
    snakemake -s up_illumina_wf.snakefile --cores 8
    ```
    The pipeline will start using 8 cores, and the results will be saved in a directory named `processed_ddmmyy` (default naming format based on the current date).  


### Project structure & key outputs

**Example project layout** after cloning the repo and executing the workflow:

```
updated_illumina_workflow/
├── raw_data/                            # SET-UP BY USER [IMPORTANT]
│   ├── run123/                          # Each run has its OWN FOLDER: runXYZ
│   │    ├── sampleA_R1_001.fastq.gz     # paired-end reads for each sample
│   │    ├── sampleA_R2_001.fastq.gz     # strictly follow naming convention
│   │    ├── sampleB_R1_001.fastq.gz     # (can have) several samples per run
│   │    └── ...
│   └── run456/                          # run names MUST be UNIQUE
│        ├── sampleC_R1_001.fastq.gz     # sample names don't have to be unique
│        └── ...                         # across the runs
├── up_illumina_wf.snakefile             # snakefile
├── environment.yaml                     # dependencies for conda installation
├── execute-and-log.sh                   # easy-to-use wrapper
├── README.md                          
└── processed_ddmmyy/                    # RESULTS: AUTO generated by snakemake
    ├── run123/
    │   └── sampleA/
    │       ├── completed_sampleA_annotation.tsv
    │       ├── sampleA_coverage.txt
    │       ├── raw/
    │       ├── dedup_qc/
    │       ├── filtered/
    │       ├── assembly/
    │       ├── mappings/
    │       ├── summary/
    │       └──  ...
    ├── run123/
    │   └── sampleC/
    │       ├──  ...
    └── ...
```

- `raw_data/` - Contains the raw FASTQ files to be processed.
- `processed_ddmmyy/` - Generated by the workflow, this directory stores processed outputs. Files are organized into subfolders for each `{run}` and `{sample}`, following the directory structure `{OUTPUT_FOLDER}/{run}/{sample}/{rule_output}/`  


### Output directories & files

1. `completed_{sample}_annotation.tsv` - Stores merged `DIAMOND` BLASTX annotation results for downstream analysis.  

2. `dedup_qc/` - Contains the **QC'ed** and **deduplicated** reads, along with `fastp` reports.  

3. `filtered/` - Holds the **human-filtered** reads (after removing host contamination.)  

4. `assembly/contigs.fasta` - Contains the contigs generated from metaSPAdes **assembly**.  

5. `{sample}_coverage.txt` - Stores the depth at each position of each contig.

6. `mappings/` -  Contains `{sample}_mappings.bam` (a binary file with filtered-reads mapped back to contigs) and `{sample}_readmap.txt` (a human-readable file showing reads that map to each contig).

7. `summary/` - `{run}_{sample}_read_processing_stats.tsv` stores read statistics across all stages of read processing; `{sample}_bam_stat.tsv` stores the total number of alignments (including multiple alignments of the same read).


## Advanced usage and configuration

The workflow can be launched by specifying the number of cores, memory, or other parameters.  
```
snakemake -s up_illumina_wf.snakefile \
    --resources mem_gb=192 \
    --cores 24 \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --latency-wait 30 \
    --max-jobs-per-second 2 \
    --max-status-checks-per-second 4
```  

- **Default output**: If no custom folder is specified by passing the `config` flag, default output directory is `processed_ddmmyy`.
- `cores`: Use up to 24 CPU cores
- `resources mem_gb`: Allocate 192 GB memory (roughly 8× the number of cores).  
- `rerun-triggers`: Force a rerun if file modification times indicate that inputs have changed.  
- `rerun-incomplete`: Rerun incomplete jobs from previous executions.
- `latency-wait`: Waits up to 30 seconds for input files (useful for network filesystems).
- **Job submission and status check rate**: Limit new job submission rate (`max-jobs-per-second`) to 2 new jobs per second; limit job status check rate (`max-status-checks-per-second`) to 4 checks per second  
- **Configuration**: Pass configuration variable flag to set a custom output folder:

  ```
  snakemake -s up_illumina_wf.snakefile \
      --config OUTPUT_FOLDER="processed_mydate" \
      --cores 16
  ```  


### Resource Allocation and Priorities  
  - Threads and Memory: Each rule in the Snakefile can dynamically allocate threads and memory.  
  - Adjusting Resources: Modify the `threads:` or `resources:` directives inside each rule for finer control.
  - Adjusting Priorities: Some rules have **priority settings** to optimize execution order. This ensures that computationally intensive steps start earlier, preventing bottlenecks and minimizing total runtime. A higher priority means the rule is executed first.
  **Example Priority Assignments**:  
    - `blastx_assembled` has **priority 2** → The most time-consuming step, so it runs first.  
    - `assemble_filtered` has **priority 1** → Runs before all other jobs except `blastx_assembled`.  

    > **Note:** By default, all rules are priority 0. Adjust priority levels as needed to optimize execution.


## Acknowledgements  

**Adapted by:** [Div Prasad](https://github.com/divprasad/) (Jul'24–Feb'25)  
**Original Workflow by:** Nathalie Worp & David Nieuwenhuijse  

This `snakemake` pipeline is an adaptation of the original work by **Nathalie Worp** and **David Nieuwenhuijse**, incorporating updates, enhancements, and additional features. Special thanks to both of them for their contributions, as their work laid the foundation for this repository.  
