#!/bin/bash

set -euo pipefail  # Enable strict error handling

# `set -euo pipefail` a commonly used safety measure in shell scripting that helps catch errors early and prevents silent failures. Here’s what each flag does:
# `-e` (Exit on Error)
#
#   Causes the script to exit immediately if any command fails (returns a non-zero exit code).
#   This helps prevent the script from continuing execution after an error, which could lead to unintended behavior.
#
# `-u` (Unset Variables Treated as Errors)
#   Treats unset variables as errors and causes the script to fail if it tries to use one.
#   This prevents issues caused by typos or missing variables.
#
# `-o` pipefail (Fail Pipeline on Error)
#
#   By default, in a pipeline (cmd1 | cmd2 | cmd3), only the exit status of the last command is considered.
#   With pipefail, if any command in the pipeline fails, the entire pipeline fails.
#
# Why Use set -euo pipefail?
#   Ensures fail-fast behavior: If a command fails, the script exits immediately instead of continuing with potential errors.
#   Catches unintended issues, like unset variables.
#   Helps debug pipelines where errors might otherwise be ignored.

# Get the current date in DDMMYY format
var_date_time=$(date +%d%m%y) # $var_date_time is now ddmmyy

# Define log file with timestamped name
log_file="logging_${var_date_time}"

# Prompt user for number of CPU cores, with a default value of 24
read -p "Enter number of cores (default 24): " cores
cores=${cores:-24}

# Prompt user for memory in GB, with a default value of 192
read -p "Enter memory in GB for resources (default 192): " mem_gb
mem_gb=${mem_gb:-192}

# Prompt user for Snakemake mode, default is dry run (-n)
read -p "Enter mode (default: dry run [n], execute [e]): " mode
mode=${mode:-n}

# Prompt user for custom output folder name, default is NO
read -p "Provide custom output folder name? (default: NO, or enter a folder name): " custom_folder_name
custom_folder_name=${custom_folder_name:-NO}

# Validate custom folder name (allow only letters, numbers, underscores, and dashes)
if [[ "$custom_folder_name" != "NO" && ! "$custom_folder_name" =~ ^[a-zA-Z0-9_-]+$ ]]; then
    echo "Error: Custom folder name contains spaces or illegal characters!"
    echo "Allowed: Letters (A-Z, a-z), Numbers (0-9), Underscores (_), and Dashes (-)."
    exit 1  # Exit script with error
fi

# Set output folder name if custom folder is provided
if [[ "$custom_folder_name" != "NO" ]]; then
    log_file="logging_${custom_folder_name}"
fi

# Print start information to the console and log file
echo "##======STARTING WORKFLOW======##" | tee -a "$log_file"
date | tee -a "$log_file"
echo "##=============================##" | tee -a "$log_file"

# Print the Snakemake execution mode
if [[ "$mode" == "e" ]]; then
    echo "Executing workflow: 'up_illumina_wf.snakefile' with ${cores} cores and ${mem_gb} GB memory" | tee -a "$log_file"
else
    echo "Performing dry run (-n): 'up_illumina_wf.snakefile' with ${cores} cores and ${mem_gb} GB memory" | tee -a "$log_file"
fi

# Define /usr/bin/time command to track execution time
usr_time="/usr/bin/time"

# Set Snakemake command options
snakemake_cmd="snakemake -s up_illumina_wf.snakefile \
    --resources mem_gb=${mem_gb} \
    --cores ${cores} \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --latency-wait 30 \
    --max-jobs-per-second 2 \
    --max-status-checks-per-second 4"

# Append -n flag if mode is dry run
if [[ "$mode" != "e" ]]; then
    snakemake_cmd+=" -n"
fi

# Set output folder name if custom folder is provided
if [[ "$custom_folder_name" != "NO" ]]; then
    snakemake_cmd+=" --config OUTPUT_FOLDER=$custom_folder_name"
fi

# Execute Snakemake
${usr_time} -o "$log_file" --append $snakemake_cmd

# Print end information to log file
echo "##======WORKFLOW COMPLETED======##" | tee -a "$log_file"
date | tee -a "$log_file"
echo -e "##=============================##\n" | tee -a "$log_file"
