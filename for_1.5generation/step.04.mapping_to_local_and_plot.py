import os

# ------------------------------------------------------------------->>>>>>>>>>
# Author: Herman ZHAO
# Email: hermanzhaozzzz@gmail.com
# Update
#     date: 2021-08-03 log: format annotation and software path
# ------------------------------------------------------------------->>>>>>>>>>

# ------------------------------------------------------------------->>>>>>>>>>
# SAMPLE INFO
# ------------------------------------------------------------------->>>>>>>>>>
SAMPLES = [
    'BAN10_A04_A04_i10',  # 2
    'BAN10_B04_B04_i10',  # 4
    'BAN10_C04_C04_i10',  # 10
    'BAN10_D04_D04_i10',  # 12
    'BAN10_E04_E04_i10',  # si
    'BAN10_F04_F04_i10'   # ba
]

READ_IDX = ["1"]
# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
THREADS = '20'

# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# polaris
# PYTHON = "/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/miniconda3/envs/snakepipes_py27/bin/python"
# abyss
# PYTHON = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq/bin/python"

dt = {
    # "ND1-on-target": "ACTCAATCCTCTGATC",
    # "ND4-on-target": "CCTGATCAAATATC",
    # "ND5.1-on-target": "ACCTTTCCTCACAGGT",
    # "ND5.3-on-target": "CTACTCATCTTCCTAATT",
    "ND6": "CCTCAGGATACTCCT"
}

# get the application path
# with os.popen("which cutadapt") as path:
#     CUTADAPT = path.read().strip()
#     print('PATH cutadapt:', CUTADAPT)
# with os.popen("which bwa") as path:
#     BWA = path.read().strip()
#     print('PATH bwa:', BWA)
# with os.popen("which bowtie") as path:
#     BOWTIE = path.read().strip()
#     BOWTIE_BUILD = BOWTIE + '-build'
#     print('PATH bowtie:', BOWTIE)
#     print('PATH bowtie-build:', BOWTIE_BUILD)
# with os.popen("which samtools") as path:
#     SAMTOOLS = path.read().strip()
#     print('PATH samtools:', SAMTOOLS)    
# with os.popen("which bedtools") as path:
#     BEDTOOLS = path.read().strip()
#     print('PATH bedtools:', BEDTOOLS)    
# with os.popen("which samclip") as path:
#     SAMCLIP = path.read().strip()
#     print('PATH samclip:', SAMCLIP)
    
# if "" in [BWA, CUTADAPT, BOWTIE, SAMTOOLS, SAMCLIP]:
#     raise ValueError("The necessary software path is missing!!!")
# else:
#     print("pass software check...\n")
CUTADAPT = "~/0.apps/micromamba/envs/snakepipes_target-seq/bin/cutadapt"
BWA = "~/0.apps/micromamba/envs/snakepipes_target-seq/bin/cutadapt/bwa"

rule all:
    input:
        expand("../../1.5generation/mapping/{sample}_R1_cutadapt.fq.gz",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.bam",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.bam.bai",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.mpileup",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.bmat",sample=SAMPLES),
        expand("../../all_plot/1.5generation/{sample}.ext50.pdf",sample=SAMPLES),
rule cutadapt:
    input:
        "../../fastq/{sample}_R1.fastq.gz"
    output:
        "../../1.5generation/mapping/{sample}_R1_cutadapt.fq.gz"
    log:
        "../../1.5generation/mapping/{sample}_cutadapt.log"
    shell:
        """
        {CUTADAPT} -j {THREADS} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -o {output} {input} > {log} 2>&1
        """
rule bwa_mapping:
    input:
        "../../1.5generation/mapping/{sample}_R1_cutadapt.fq.gz"
    output:
        temp("../../1.5generation/mapping/{sample}_bwa.sam")
    params:
        ref = '../../reference.fasta/plasmid.ref.upper.fa'
    log:
        "../../1.5generation/mapping/{sample}_bwa.log"
    shell:
        """
        bwa mem {params.ref} {input} -t {THREADS} > {output} 2>{log}
        """
rule sam_to_bam:
    input:
        sam = "../../1.5generation/mapping/{sample}_bwa.sam",
    output:
        temp("../../1.5generation/mapping/{sample}_bwa.bam")
    shell:
        """
        samtools sort -n {input.sam} > {output}
        """
rule samtools_sort_by_position:
    input:
        "../../1.5generation/mapping/{sample}_bwa.bam"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam"
    shell:
        "samtools sort -O BAM -o {output} -T {output}.temp -@ {THREADS} {input}"
rule samtools_index:
    input:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam.bai"
    shell:
        "samtools index -@ 4 {input} {output}"
rule spike_in_mpileup:
    input:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam",
        "../../1.5generation/mapping/{sample}_bwa_sort.bam.bai"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.mpileup"
    params:
        ref = '../../reference.fasta/plasmid.ref.upper.fa'
    shell:
        "samtools mpileup {input[0]} --reference {params.ref} --max-depth 10000000 > {output} " 
rule parse_mpileup:
    input: 
        "../../1.5generation/mapping/{sample}_bwa_sort.mpileup"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.bmat"
    shell:
        "bioat bam mpileup_to_table {input} -o {output} --mutation_number_threshold 0"
rule bmat_plot:
    input: 
        "../../1.5generation/mapping/{sample}_bwa_sort.bmat"
    output:
        "../../all_plot/1.5generation/{sample}.ext50.pdf"
    params:
        sgRNA_seq = "CAACAAAGCCACGTTGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGAACAATAAAACTGTCTGCTTACATAAACAGTAATACAAGGGGTGTTATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCA"
    shell:
        """
        bioat target_seq region_heatmap --input_table {input} --output_fig {output} --region_extend_length 50 --target_seq {params.sgRNA_seq}
        """