# adapted by Div Prasad, written between Jul'24 - Feb'25
# updated the illumina workflow developed by Nathalie Worp and David Nieuwenhuijse

import multiprocessing
import datetime

# Define the output folder dynamically based on current day and month.
OUTPUT_FOLDER = config.get("OUTPUT_FOLDER", f"processed_{datetime.datetime.now().strftime('%d%m%y')}")

# Get the total number of cores specified when launching snakemake
total_cores = workflow.cores

# Capture the RUN and SAMPLE wildcards from the file paths
RUNS, SAMPLES = glob_wildcards("raw_data/{run}/{sample}_R1_001.fastq.gz")

rule all:
    input:
        expand(f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/assembly/contigs.fasta", zip, run=RUNS, sample=SAMPLES),
        expand(f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/summary/{{sample}}_bam_stat.tsv", zip, run=RUNS, sample=SAMPLES),
        expand(f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/summary/{{run}}_{{sample}}_read_processing_stats_clean.tsv", zip, run=RUNS, sample=SAMPLES),
        expand(f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/renamed_completed_{{sample}}_annotation.tsv", zip, run=RUNS, sample=SAMPLES),
        expand(f"{OUTPUT_FOLDER}/{{run}}/annotations/renamed/renamed_completed_{{sample}}_annotation.tsv", zip, run=RUNS, sample=SAMPLES)

rule hardlink_raw:
    input:
        R1 = "raw_data/{run}/{sample}_R1_001.fastq.gz",
        R2 = "raw_data/{run}/{sample}_R2_001.fastq.gz"
    output:
        R1 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/raw/{{run}}_{{sample}}_R1.fastq.gz",
        R2 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/raw/{{run}}_{{sample}}_R2.fastq.gz"
    shell:
        """
        set -euo pipefail

        # create hard links so to help back trace which input files have been processed
        ln $(realpath {input.R1}) {output.R1}
        ln $(realpath {input.R2}) {output.R2}
        """

rule QC_after_dedup:
    input:
        R1 = "raw_data/{run}/{sample}_R1_001.fastq.gz",
        R2 = "raw_data/{run}/{sample}_R2_001.fastq.gz"
    output:
        R1 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_R1.fastq.gz",
        R2 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_R2.fastq.gz",
        S  = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_S.fastq.gz",
        failed = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_fail.fastq.gz"
    threads: max(4, round(total_cores / 4))
    shell:
        """
        set -euo pipefail

        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} \
        --unpaired1 {output.S} --unpaired2 {output.S} --failed_out {output.failed} \
        -w {threads} --detect_adapter_for_pe --dedup \
        -l 50 -Q -y -x -c \
        --cut_right \
        --cut_right_window_size 6 \
        --cut_right_mean_quality 20 \
        -j {OUTPUT_FOLDER}/{wildcards.run}/{wildcards.sample}/dedup_qc/{wildcards.run}_{wildcards.sample}_qc_report.json \
        -h {OUTPUT_FOLDER}/{wildcards.run}/{wildcards.sample}/dedup_qc/{wildcards.run}_{wildcards.sample}_qc_report.html
        """

# -aY . Y is use softclipping -a is output all alignments
# the output of bwa-mem2 (is exactly same as bwa) mem is a sam file that is piped directly into samtools in order to not save the memory consuming sam files. The samtools fastq option -c means the level of compression when writing .gz fastq files. The output bam/sam file contains information about mapped and unmapped reads.
# -f 4 : then you obtain the unmapped (which are not human reads) reads and keep these.
#-s option: If a singleton file is specified using the -s option then only paired sequences will be output for categories 1 and 2; paired meaning that for a given QNAME there are sequences for both category 1 and 2. If there is a sequence for only one of categories 1 or 2 then it will be diverted into the specified singletons file. This can be used to prepare fastq files for programs that cannot handle a mixture of paired and singleton reads.

rule filter_human:
    input:
        R1 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_R1.fastq.gz",
        R2 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_R2.fastq.gz",
        S  = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_S.fastq.gz"
    output:
        R1 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R1.fastq.gz",
        R2 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R2.fastq.gz",
        S  = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_S.fastq.gz"
    threads: max(4, round(total_cores / 4))
    shell:
        """
        set -euo pipefail

        bwa-mem2 mem -aY -t {threads} /mnt/scratch2/DB/HG38/bwa-mem2/GCF_000001405.26_GRCh38_genomic.fna {input.R1} {input.R2} | \
        samtools fastq -c 3 -f 4 -1 {output.R1} -2 {output.R2} -s /dev/null -
        bwa-mem2 mem -aY -t {threads} /mnt/scratch2/DB/HG38/bwa-mem2/GCF_000001405.26_GRCh38_genomic.fna {input.S} | \
        samtools fastq -c 3 -f 4 - | gzip > {output.S}
        """

#de novo assembly of contigs
#file with forward paired-end reads + file with reverse paired end reads
# -s file with unpaired reads

rule assemble_filtered:
    priority: 1
    input:
        R1 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R1.fastq.gz",
        R2 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R2.fastq.gz",
        S  = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_S.fastq.gz"
    output:
        assembled_contigs    = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/assembly/contigs.fasta",
        renamed_contigs      = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/assembly/{{run}}_{{sample}}_contigs.fasta"
    params:
        spades_folder = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/assembly"
    threads: max(16, round(total_cores / 4))
    shell:
        """
        set -euo pipefail

        metaspades.py -t {threads} \
        --tmp-dir /dev/shm \
        -o {params.spades_folder} \
        -k 21,33,45,55,67,77,99,113,127 \
        -1 {input.R1} \
        -2 {input.R2} \
        -s {input.S}

        sed "s/^>/>{wildcards.run}_{wildcards.sample}|/" {output.assembled_contigs} > {output.renamed_contigs}
        """
        # Append the filename to the beginning of each contig

# tantan: identify low-complexity regions in the sequences.
# Mask these regions (make them ignored in alignment) while leaving the regions visible in the sequence for reference.

rule blastx_assembled:
    priority: 2
    input:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/assembly/contigs.fasta"
    output:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/{{sample}}_annotation.tsv"
    threads: max(16, round(total_cores / 3))
    shell:
        """
        set -euo pipefail

        diamond blastx \
        -q {input} \
        -o {output} \
        -f 6 qseqid qlen sseqid slen approx_pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle \
        -d /mnt/scratch2/DB/dmnd/dmnd_nr_db_2024-07-01/2024-07-01_nr_db.dmnd \
        --tmpdir /dev/shm --evalue 10-5 --max-target-seqs 50 --soft-masking tantan \
        --threads {threads} -c1 -b15.0 --fast
        """

# map reads back to the assembled contigs; needed to calculate the depth per position, mean coverage per contig etc
rule map_reads_to_contigs:
    input:
        R1      = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R1.fastq.gz",
        R2      = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R2.fastq.gz",
        S       = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_S.fastq.gz",
        contigs = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/assembly/contigs.fasta"
    params:
        tmp_paired   = temp("{sample}_paired.tmp"),
        tmp_singlets = temp("{sample}_singlets.tmp")
    output:
        bam = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/mappings/{{sample}}_mappings.bam"
    threads: max(4, round(total_cores / 4))
    shell:
        """
        set -euo pipefail

        bwa-mem2 index {input.contigs}
        bwa-mem2 mem -Y -t {threads} {input.contigs} {input.R1} {input.R2} | samtools view -b - | samtools sort -o {params.tmp_paired}
        bwa-mem2 mem -Y -t {threads} {input.contigs} {input.S} | samtools view -b - | samtools sort -o {params.tmp_singlets}
        samtools merge {output.bam} {params.tmp_paired} {params.tmp_singlets}
        """

# -b: output in de bam format, -F2052 Do not output alignments with any bits set in INT present in the FLAG field. INT = 2052 means that alignments that are not mapped are not outputted.
# samtools depth: Computes the depth at each position or region. -a outputs all positions (including those with zero depth) -d
# -d 0: read at most 0 reads per input file. This means figures greater than 0 may be reported in the output. Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit.

rule create_coverage_file:
    input:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/mappings/{{sample}}_mappings.bam"
    output:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/{{sample}}_coverage.txt"
    threads: 1
    shell:
        """
        set -euo pipefail
        samtools view -bF2052 {input} | samtools depth -a -d 0 - > {output}
        """

#, do not output unmapped alignments and print the first and third field to output
rule create_readmap_file:
    input:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/mappings/{{sample}}_mappings.bam"
    output:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/mappings/{{sample}}_readmap.txt"
    shell:
        """
        set -euo pipefail
        samtools view -F2052 {input} | cut -f1,3 > {output}
        """


rule stat_mapped_files:
    input:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/mappings/{{sample}}_mappings.bam"
    output:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/summary/{{sample}}_bam_stat.tsv"
    shell:
        """
        set -euo pipefail
        samtools view -c -F2052 {input} > {output}
        """

# combine all the classification and coverage results to create summary tables
rule merge_results:
    input:
        annotation = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/{{sample}}_annotation.tsv",
        coverage   = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/{{sample}}_coverage.txt",
        readmap    = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/mappings/{{sample}}_readmap.txt"
    output:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/completed_{{sample}}_annotation.tsv"
    script:
        "./merge_results.R"


# create links to all the "completed*" classification tables per sample, to a central location per run
# Skip the first line, then prepend "{run}_{sample}|" to each remaining line
# Remove temp files associated with the sample
rule store_completed_annotation_files:
    input:
        f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/completed_{{sample}}_annotation.tsv"
    output:
        renamed_completed    = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/renamed_completed_{{sample}}_annotation.tsv",
        ln_completed         = f"{OUTPUT_FOLDER}/{{run}}/annotations/completed_{{sample}}_annotation.tsv",
        ln_renamed_completed = f"{OUTPUT_FOLDER}/{{run}}/annotations/renamed/renamed_completed_{{sample}}_annotation.tsv"
    shell:
        """
        set -euo pipefail

        sed "1!s/^/{wildcards.run}_{wildcards.sample}|/" {input} > {output.renamed_completed}

        cp {input} {output.ln_completed}
        cp {output.renamed_completed} {output.ln_renamed_completed}

        rm -f "{wildcards.run}_{wildcards.sample}"*
        rm -f "{wildcards.sample}"*

        """

        # # For easier access, create hard links instead of copying files
        # ln $(realpath {input}) {output.ln_completed}
        # ln $(realpath {output.renamed_completed}) {output.ln_renamed_completed}


rule seqkit_stat:
    input:
        raw_R1          = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/raw/{{run}}_{{sample}}_R1.fastq.gz",
        raw_R2          = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/raw/{{run}}_{{sample}}_R2.fastq.gz",
        qc_dedup_R1     = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_R1.fastq.gz",
        qc_dedup_R2     = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_R2.fastq.gz",
        qc_dedup_S      = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/dedup_qc/{{run}}_{{sample}}_S.fastq.gz",
        fil_qc_dedup_R1 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R1.fastq.gz",
        fil_qc_dedup_R2 = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_R2.fastq.gz",
        fil_qc_dedup_S  = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/filtered/{{run}}_{{sample}}_S.fastq.gz"
    output:
        combined_stats       = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/summary/{{run}}_{{sample}}_read_processing_stats.tsv",
        clean_combined_stats = f"{OUTPUT_FOLDER}/{{run}}/{{sample}}/summary/{{run}}_{{sample}}_read_processing_stats_clean.tsv"
    shell:
        """
        set -euo pipefail

        seqkit stat {input.raw_R1} {input.raw_R2} \
        {input.qc_dedup_R1} {input.qc_dedup_R2} {input.qc_dedup_S} \
        {input.fil_qc_dedup_R1} {input.fil_qc_dedup_R2} \
        {input.fil_qc_dedup_S} |\
        awk 'FNR==1 && NR!=1 {{next}}{{print}}' \
        > {output.combined_stats}

        sed 's/,//g' {output.combined_stats} | sed 's/  */ /g' \
        > {output.clean_combined_stats}
        """
        # clean_combined_stats essentially removes left trailing white spaces,
        # replaces tabs with single space (" ")
        # and removes the commas for various numbers
        # this helps to easily handle the `clean_combined_stats` for downstream analysis
