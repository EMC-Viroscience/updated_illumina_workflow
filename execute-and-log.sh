#!/bin/bash

set -euo pipefail  # Enable strict error handling

"""
`set -euo pipefail` a commonly used safety measure in shell scripting that helps catch errors early and prevents silent failures. Here’s what each flag does:

`-e` (Exit on Error)

  Causes the script to exit immediately if any command fails (returns a non-zero exit code).
  This helps prevent the script from continuing execution after an error, which could lead to unintended behavior.

`-u` (Unset Variables Treated as Errors)
  Treats unset variables as errors and causes the script to fail if it tries to use one.
  This prevents issues caused by typos or missing variables.

`-o` pipefail (Fail Pipeline on Error)

  By default, in a pipeline (cmd1 | cmd2 | cmd3), only the exit status of the last command is considered.
  With pipefail, if any command in the pipeline fails, the entire pipeline fails.

*Why Use set -euo pipefail?*
  Ensures fail-fast behavior: If a command fails, the script exits immediately instead of continuing with potential errors.
  Catches unintended issues, like unset variables.
  Helps debug pipelines where errors might otherwise be ignored.
"""


# Get the current date in DDMMYY format
var_date_time=$(date +%d%m%y)
# $var_date_time is now ddmmyy

# Define log file with timestamped name
log_file="logging_${var_date_time}"

# Remove any previous logs (optional, uncomment if needed)
# rm -rf logging_*

# Prompt user for number of CPU cores, with a default value of 16
read -p "Enter number of cores (default 16): " cores
cores=${cores:-24}

# Prompt user for memory in GB, with a default value of 128
read -p "Enter memory in GB for resources (default 128): " mem_gb
mem_gb=${mem_gb:-192}

# Print start information to the console and log file
echo "##======STARTING WORKFLOW======##" | tee -a "$log_file"
date | tee -a "$log_file"
echo "##=============================##" | tee -a "$log_file"

# Inform the user about the Snakemake execution
echo "Running 'updated_illumina_workflow.smk' with ${cores} cores and ${mem_gb} GB memory" | tee -a "$log_file"

# Define /usr/bin/time command to track execution time
usr_time="/usr/bin/time"

# Execute the Snakemake workflow with resource limits
# -n flag is for dry run (doesn't execute the workflow, but to check if errors in DAG)
# once dry run jobs make sense, remove -n flag

${usr_time} -o "$log_file" --append \
    snakemake -n -s updated_illumina_workflow.smk \  # Specify the workflow file
    --config OUTPUT_FOLDER="processed_${var_date_time}" \  # Set the output directory
    # by default snakemake will set output directory as "processed_ddmmyy"
    --resources mem_gb=${mem_gb} \  # Allocate the specified GB of memory
    --cores ${cores} \  # Use the specified number of CPU cores
    --rerun-triggers mtime \  # Force rerun if input files are modified (based on modification time)
    --rerun-incomplete \  # Retry any incomplete jobs from previous runs
    --latency-wait 30 \  # Wait up to 30 seconds for input files (useful for network filesystems)
    --max-jobs-per-second 2 \  # Limit job submission rate to 2 per second
    --max-status-checks-per-second 4  # Limit status checks to 4 per second

# Print end information to log file
echo "##======WORKFLOW COMPLETED======##" | tee -a "$log_file"
date | tee -a "$log_file"
echo "##=============================##" | tee -a "$log_file"
