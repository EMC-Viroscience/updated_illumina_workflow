#!/bin/bash

# Prompt user for cores and memory resources with defaults
read -p "Enter number of cores (default 16): " cores
cores=${cores:-16}

read -p "Enter memory in GB for resources (default 128): " mem_gb
mem_gb=${mem_gb:-128}

# Print start information to the console
echo "##------START------##"
date
echo ""

# Note: The --rerun-triggers mtime option forces Snakemake to rerun jobs
# if any input file has a modification time that is more recent than its output.
echo "Running 'updated_illumina_workflow_241003.smk' with ${cores} cores and ${mem_gb} GB memory"

# Execute the Snakemake workflow.
# Standard output is discarded and only errors are appended to log.txt.
# /usr/bin/time appends timing information to log.txt.
usr_time="/usr/bin/time"
${usr_time} -o log.txt --append snakemake -s updated_illumina_workflow_241003.smk \
    --rerun-triggers mtime \  # Forces rerun if input files are updated
    --rerun-incomplete \
    --cores ${cores} --resources mem_gb=${mem_gb} \
    --latency-wait 30 --max-jobs-per-second 2 --max-status-checks-per-second 4 \
    > /dev/null 2>> log.txt

date
echo "##------END------##"
