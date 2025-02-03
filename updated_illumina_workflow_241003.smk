# adapted by Divyae Kishore Prasad, written between July-October 2024
# updated the "illumina workflow" developed by Nathalie Worp and David Nieuwenhuisje

import multiprocessing

# Get the total number of cores specified when launching snakemake
total_cores = workflow.cores

# Capture the RUN and SAMPLE wildcards from the file paths
RUNS, SAMPLES = glob_wildcards("raw_data/{run}/{sample}_R1_001.fastq.gz")

rule all:
    input:
        expand("processed_0209/{run}/{sample}/assembly/contigs.fasta", zip, run=RUNS, sample=SAMPLES),
        expand("processed_0209/{run}/{sample}/summary/{sample}_bam_stat.tsv", zip, run=RUNS, sample=SAMPLES),
        expand("processed_0209/{run}/{sample}/summary/{sample}_readstats_raw.tsv", zip, run=RUNS, sample=SAMPLES),
        expand("processed_0209/{run}/{sample}/summary/{sample}_readstats_dedup_qc.tsv", zip, run=RUNS, sample=SAMPLES),
        expand("processed_0209/{run}/{sample}/summary/{sample}_readstats_filtered.tsv", zip, run=RUNS, sample=SAMPLES),
        expand("processed_0209/{run}/{sample}/summary/{run}_{sample}_read_processing_stats_clean.tsv", zip, run=RUNS, sample=SAMPLES),
        expand("processed_0209/{run}/{sample}/renamed_completed_{sample}_annotation.tsv", zip, run=RUNS, sample=SAMPLES),
        expand("processed_0209/{run}/annotations/renamed/renamed_completed_{sample}_annotation.tsv", zip, run=RUNS, sample=SAMPLES)

rule softlink_raw:
    input:
        R1 = "raw_data/{run}/{sample}_R1_001.fastq.gz",
        R2 = "raw_data/{run}/{sample}_R2_001.fastq.gz"
    output:
        R1 = "processed_0209/{run}/{sample}/raw/{run}_{sample}_R1.fastq.gz",
        R2 = "processed_0209/{run}/{sample}/raw/{run}_{sample}_R2.fastq.gz"
    shell:
        """
        set -euo pipefail

        ln -s $(realpath {input.R1}) {output.R1}
        ln -s $(realpath {input.R2}) {output.R2}
        """

rule QC_after_dedup:
    input:
        R1 = "raw_data/{run}/{sample}_R1_001.fastq.gz",
        R2 = "raw_data/{run}/{sample}_R2_001.fastq.gz"
    output:
        R1 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R1.fastq.gz",
        R2 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R2.fastq.gz",
        S = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_S.fastq.gz",
        failed = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_fail.fastq.gz"
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
        -j processed_0209/{wildcards.run}/{wildcards.sample}/dedup_qc/{wildcards.run}_{wildcards.sample}_qc_report.json \
        -h processed_0209/{wildcards.run}/{wildcards.sample}/dedup_qc/{wildcards.run}_{wildcards.sample}_qc_report.html
        """

# changed FROM -> -l 30 -Q -y --cut_right_window_size 5 --cut_right_mean_quality 25 \
# for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it.
# so file S contains reads. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])

rule filter_human:
    input:
        R1 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R1.fastq.gz",
        R2 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R2.fastq.gz",
        S = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_S.fastq.gz"
    output:
        R1 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R1.fastq.gz",
        R2 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R2.fastq.gz",
        S = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_S.fastq.gz"
    threads: max(4, round(total_cores / 4))
    shell:
        """
        set -euo pipefail

        bwa-mem2 mem -aY -t {threads} /mnt/scratch2/DB/HG38/bwa-mem2/GCF_000001405.26_GRCh38_genomic.fna {input.R1} {input.R2} | \
        samtools fastq -c 3 -f 4 -1 {output.R1} -2 {output.R2} -s /dev/null -
        bwa-mem2 mem -aY -t {threads} /mnt/scratch2/DB/HG38/bwa-mem2/GCF_000001405.26_GRCh38_genomic.fna {input.S} | \
        samtools fastq -c 3 -f 4 - | gzip > {output.S}
        """
# -aY . Y is use softclipping -a is output all alignments
# the output of bwa-mem2 (is exactly same as bwa) mem is a sam file that is piped directly into samtools in order to not save the memory consuming sam files. The samtools fastq option -c means the level of compression when writing .gz fastq files. The output bam/sam file contains information about mapped and unmapped reads.
# -f 4 : then you obtain the unmapped (which are not human reads) reads and keep these.
#-s option: If a singleton file is specified using the -s option then only paired sequences will be output for categories 1 and 2; paired meaning that for a given QNAME there are sequences for both category 1 and 2. If there is a sequence for only one of categories 1 or 2 then it will be diverted into the specified singletons file. This can be used to prepare fastq files for programs that cannot handle a mixture of paired and singleton reads.

rule assemble_filtered:
    input:
        R1 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R1.fastq.gz",
        R2 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R2.fastq.gz",
        S = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_S.fastq.gz"
    output:
        assembled_contigs = "processed_0209/{run}/{sample}/assembly/contigs.fasta",
        renamed_contigs = "processed_0209/{run}/{sample}/assembly/{run}_{sample}_contigs.fasta",
        fil_renamed_contigs = "processed_0209/{run}/{sample}/assembly/fil_{run}_{sample}_contigs.fasta"
    params:
        spades_folder = "processed_0209/{run}/{sample}/assembly"
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

        # Append the filename to the beginning of each contig
        sed "s/^>/>{wildcards.run}_{wildcards.sample}|/" {output.assembled_contigs} > {output.renamed_contigs}

        python /home/r061231/scripts/upper_FastaMLtoSL.py {output.renamed_contigs}
        mv {output.renamed_contigs}_SL.fa {output.renamed_contigs}

        seqkit seq -g -m 250 {output.renamed_contigs} > {output.fil_renamed_contigs}
        python /home/r061231/scripts/upper_FastaMLtoSL.py {output.fil_renamed_contigs}
        mv {output.fil_renamed_contigs}_SL.fa {output.fil_renamed_contigs}
        """
#de novo assembly of contigs
#file with forward paired-end reads + file with reverse paired end reads
# -s file with unpaired reads

rule blastx_assembled:
    input:
        "processed_0209/{run}/{sample}/assembly/contigs.fasta"
    output:
        "processed_0209/{run}/{sample}/{sample}_annotation.tsv"
    threads: max(16, round(total_cores / 3))
    resources:
        mem_gb=round(360//3)  # Specify the amount of memory in GB, adjust accordingly
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
# tantan: identify low-complexity regions in the sequences.
# Mask these regions (make them ignored in alignment) while leaving the regions visible in the sequence for reference.
# --taxonmap /mnt/viro0002/workgroups_projects/Bioinformatics/DB/nr_2024-06-18/prot.accession2taxid.FULL.gz \

# map reads back to the assembled contigs; needed to calculate the depth per position, mean coverage per contig etc
rule map_reads_to_contigs:
    input:
        R1 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R1.fastq.gz",
        R2 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R2.fastq.gz",
        S = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_S.fastq.gz",
        contigs = "processed_0209/{run}/{sample}/assembly/contigs.fasta"
    params:
        tmp_paired = temp("{sample}_paired.tmp"),
        tmp_singlets = temp("{sample}_singlets.tmp")
    output:
        bam = "processed_0209/{run}/{sample}/mappings/{sample}_mappings.bam"
    threads: max(4, round(total_cores / 4))
    shell:
        """
        set -euo pipefail

        bwa-mem2 index {input.contigs}
        bwa-mem2 mem -Y -t {threads} {input.contigs} {input.R1} {input.R2} | samtools view -b - | samtools sort -o {params.tmp_paired}
        bwa-mem2 mem -Y -t {threads} {input.contigs} {input.S} | samtools view -b - | samtools sort -o {params.tmp_singlets}
        samtools merge {output.bam} {params.tmp_paired} {params.tmp_singlets}
        """

# stores the depth at each contig position, for all contigs
rule create_coverage_file:
    input:
        "processed_0209/{run}/{sample}/mappings/{sample}_mappings.bam"
    output:
        "processed_0209/{run}/{sample}/{sample}_coverage.txt"
    threads: 1
    shell:
        """
        set -euo pipefail
        samtools view -bF2052 {input} | samtools depth -a -d 0 - > {output}
        """
# -b: output in de bam format, -F2025 Do no2 output alignments with any bits set in INT present in the FLAG field. INT = 2025 and means that alignments that are not mapped are not outputted.
#samtools depth: Computes the depth at each position or region. -a Output all positions (including those with zero depth) -d
# -d 0:  read at most 0 (INT) reads per input file. This means figures greater than INT may be reported in the output.Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit.

rule create_readmap_file:
    input:
        "processed_0209/{run}/{sample}/mappings/{sample}_mappings.bam"
    output:
        "processed_0209/{run}/{sample}/mappings/{sample}_readmap.txt"
    shell:
        """
        set -euo pipefail
        samtools view -F2052 {input} | cut -f1,3 > {output}
        """
#, do not output unmapped alignments and print the first and third field to output

rule stat_mapped_files:
    input:
        "processed_0209/{run}/{sample}/mappings/{sample}_mappings.bam"

    output:
        "processed_0209/{run}/{sample}/summary/{sample}_bam_stat.tsv"
    shell:
        """
        set -euo pipefail
        samtools view -c -F2052 {input} > {output}
        """

# combine all the classification and coverage results to create summary tables
rule merge_results:
    input:
        annotation = "processed_0209/{run}/{sample}/{sample}_annotation.tsv",
        coverage = "processed_0209/{run}/{sample}/{sample}_coverage.txt",
        readmap = "processed_0209/{run}/{sample}/mappings/{sample}_readmap.txt"
    output:
        "processed_0209/{run}/{sample}/completed_{sample}_annotation.tsv"
    script:
        "/mnt/scratch2/ww_virome_capture_longitudinal/merge_results.R"

# create links to all the "completed*" classification tables per sample, to a central location per run
rule store_completed_annotation_files:
    input:
        "processed_0209/{run}/{sample}/completed_{sample}_annotation.tsv"
    output:
        renamed_completed = "processed_0209/{run}/{sample}/renamed_completed_{sample}_annotation.tsv",
        ln_completed = "processed_0209/{run}/annotations/completed_{sample}_annotation.tsv",
        ln_renamed_completed = "processed_0209/{run}/annotations/renamed/renamed_completed_{sample}_annotation.tsv"
    shell:
        """
        set -euo pipefail

        # Extract the sample number before '_S' and prepend 'bc'
        # sample_number=$(echo {wildcards.sample} | awk -F'_S' '{{print "bc"$1}}')

        # ```sed "1 ! s/^/``` skip the first line, and then insert at the beginning of each line
        sed "1 ! s/^/{wildcards.run}_{wildcards.sample}|/" {input} > {output.renamed_completed}

        # # Copy files for easier access
        ln -s $(realpath {input}) {output.ln_completed}
        ln -s $(realpath {output.renamed_completed}) {output.ln_renamed_completed}

        # # For easier access, create symbolic links instead of copying files
        # ln -s $(realpath {input}) {output.ln_completed}
        # ln -s $(realpath {output.renamed_completed}) {output.ln_renamed_completed}
        # # Remove temp files associated with the sample
        # rm -f {wildcards.run}_{wildcards.sample}*
        # rm -f {wildcards.sample}*

        """

rule raw_stats:
    input:
        r1 = "processed_0209/{run}/{sample}/raw/{run}_{sample}_R1.fastq.gz",
        r2 = "processed_0209/{run}/{sample}/raw/{run}_{sample}_R2.fastq.gz"
    output:
        "processed_0209/{run}/{sample}/summary/{sample}_readstats_raw.tsv"
    shell:
        """
        set -euo pipefail
        seqkit stats {input.r1} {input.r2} | awk 'FNR==1 && NR!=1{{next}}{{print}}' > {output}
        """
        # ``` awk 'FNR==1 && NR!=1{{next}}{{print}}' ```
        # command skips the first line of all files after the first file when concatenating multiple files, effectively removing repeated headers.

rule dedup_qc_stats:
    input:
        r1 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R1.fastq.gz",
        r2 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R2.fastq.gz"
    output:
        "processed_0209/{run}/{sample}/summary/{sample}_readstats_dedup_qc.tsv"
    shell:
        """
        set -euo pipefail
        seqkit stats {input.r1} {input.r2} | awk 'FNR==1 && NR!=1{{next}}{{print}}' > {output}
        """

rule filtered_stats:
    input:
        r1 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R1.fastq.gz",
        r2 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R2.fastq.gz"
    output:
        "processed_0209/{run}/{sample}/summary/{sample}_readstats_filtered.tsv"
    shell:
        """
        set -euo pipefail
        seqkit stats {input.r1} {input.r2} | awk 'FNR==1 && NR!=1{{next}}{{print}}' > {output}
        """

rule seqkit_stat:
    input:
        raw_R1 = "processed_0209/{run}/{sample}/raw/{run}_{sample}_R1.fastq.gz",
        raw_R2 = "processed_0209/{run}/{sample}/raw/{run}_{sample}_R2.fastq.gz",
        qc_dedup_R1 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R1.fastq.gz",
        qc_dedup_R2 = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_R2.fastq.gz",
        qc_dedup_S = "processed_0209/{run}/{sample}/dedup_qc/{run}_{sample}_S.fastq.gz",
        fil_qc_dedup_R1 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R1.fastq.gz",
        fil_qc_dedup_R2 = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_R2.fastq.gz",
        fil_qc_dedup_S = "processed_0209/{run}/{sample}/filtered/{run}_{sample}_S.fastq.gz"
    output:
        combined_stats = "processed_0209/{run}/{sample}/summary/{run}_{sample}_read_processing_stats.tsv",
        clean_combined_stats = "processed_0209/{run}/{sample}/summary/{run}_{sample}_read_processing_stats_clean.tsv"
    shell:
        """
        set -euo pipefail

        seqkit stat {input.raw_R1} {input.raw_R2} \
        {input.qc_dedup_R1} {input.qc_dedup_R2} {input.qc_dedup_S} \
        {input.fil_qc_dedup_R1} {input.fil_qc_dedup_R2} \
        {input.fil_qc_dedup_S} |\
        awk 'FNR==1 && NR!=1 {{next}}{{print}}' \
        > {output.combined_stats}

        seqkit stat {input.raw_R1} {input.raw_R2} \
        {input.qc_dedup_R1} {input.qc_dedup_R2} {input.qc_dedup_S} \
        {input.fil_qc_dedup_R1} {input.fil_qc_dedup_R2} \
        {input.fil_qc_dedup_S} |\
        awk 'FNR==1 && NR!=1 {{next}}{{print}}' |\
        sed 's/,//g' | sed 's/  */ /g' \
        > {output.clean_combined_stats}
        """
