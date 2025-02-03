echo "\n##------START------##" >> log.txt
date >> log.txt
echo "" >> log.txt

echo "running \`updated_illumina_workflow_241003.smk\` with 48 cores" >> log.txt;
/usr/bin/time -o log.txt --append snakemake -s updated_illumina_workflow_241003.smk \
--rerun-triggers mtime --rerun-incomplete \
--cores 48 --resources mem_gb=480 \
--latency-wait 30 --max-jobs-per-second 2 --max-status-checks-per-second 4;
date >> log.txt

# --rerun-triggers mtime
# --latency-wait: Controls how long Snakemake waits for output files to appear after a job completes (helps with network or filesystem delays).
# --max-jobs-per-second: Limits the rate at which Snakemake submits new jobs (prevents scheduler overload).
# --max-status-checks-per-second: Limits how often Snakemake checks the status of running jobs (reduces load on the job scheduler).

echo "\n##------END------##" >> log.txt
