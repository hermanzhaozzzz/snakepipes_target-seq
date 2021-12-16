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


CUTOFF = ["0"]

LIBS = ['ND6']


SAMPLES = ['DddA-NES-ND6-rep1',
 'DddA-NES-ND6-rep2',
 'DddIA-1-0.25-ND6-rep1',
 'DddIA-1-0.25-ND6-rep2',
 'DddIA-1-0.5-ND6-rep1',
 'DddIA-1-0.5-ND6-rep2',
 'DddIA-1-1-ND6-rep1',
 'DddIA-1-1-ND6-rep2',
 'DddIA-1-1.2-ND6-rep1',
 'DddIA-1-1.2-ND6-rep2',
 'DddIA-1-1.5-ND6-rep1',
 'DddIA-1-1.5-ND6-rep2',
 'G1309A-ND6-rep1',
 'G1309A-ND6-rep2',
 'N1308A-ND6-rep1',
 'N1308A-ND6-rep2',
 'N1367A-ND6-rep1',
 'N1367A-ND6-rep2',
 'N1368A-ND6-rep1',
 'Q1310A-ND6-rep1',
 'Q1310A-ND6-rep2',
 'TALE-NES-ND6-rep1',
 'TALE-NES-ND6-rep2',
 'UGI-NES-ND6-rep1',
 'UGI-NES-ND6-rep2',
 'WT-ND6-rep1',
 'WT-ND6-rep2',
 'untreated-rep1',
 'untreated-rep2']

READ_IDX = ["1"]
# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
THREADS = '20'

# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
PRIMER_INFO = "./primer_table/primer.txt"
GENOME_HG38 = "/lustre1/chengqiyi_pkuhpc/zhaohn/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# polaris
PYTHON = "/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/miniconda3/envs/snakepipes_py27/bin/python"
# abyss
# PYTHON = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq/bin/python"

# get the application path
with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
    print('PATH cutadapt:', CUTADAPT)
with os.popen("which bwa") as path:
    BWA = path.read().strip()
    print('PATH bwa:', BWA)
with os.popen("which bowtie") as path:
    BOWTIE = path.read().strip()
    BOWTIE_BUILD = BOWTIE + '-build'
    print('PATH bowtie:', BOWTIE)
    print('PATH bowtie-build:', BOWTIE_BUILD)
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
    print('PATH samtools:', SAMTOOLS)    
with os.popen("which bedtools") as path:
    BEDTOOLS = path.read().strip()
    print('PATH bedtools:', BEDTOOLS)    
with os.popen("which samclip") as path:
    SAMCLIP = path.read().strip()
    print('PATH samclip:', SAMCLIP)
    
if "" in [BWA, CUTADAPT, BOWTIE, SAMTOOLS, SAMCLIP]:
    raise ValueError("The necessary software path is missing!!!")
else:
    print("pass software check...\n")
rule all:
    input:
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        # expand('../reference.fasta/{sample}.ref.upper.fa',sample=SAMPLES),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../all_plot/cutoff_{cutoff}.ext50/TargetSeq-{lib}_{sample}_cutoff_{cutoff}_indel.ext50.pdf",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        # expand('../reference.fasta/{sample}.sgRNA.upper.fa.seq',sample=SAMPLES)
rule cutadapt:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz"
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_cutadapt.log"
    shell:
        """
        {CUTADAPT} -j {THREADS} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -o {output} {input} > {log} 2>&1
        """
rule bwa_mapping:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam"
    params:
        ref = '../reference.fasta/WT.ref.upper.fa'
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.log"
    shell:
        """
        {BWA} mem {params.ref} {input} -t {THREADS} -M > {output} 2>{log}
        """
rule sam_to_bam:
    input:
        sam = "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam",
        ref = '../reference.fasta/WT.ref.upper.fa'
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    shell:
        """
        {SAMTOOLS} view -hb {input.sam}> {output}
        """
rule samtools_sort_by_position:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    shell:
        "{SAMTOOLS} sort -O BAM -o {output} -T {output}.temp -@ {THREADS} -m 2G {input}"
rule samtools_index:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai"
    shell:
        "{SAMTOOLS} index -@ {THREADS} {input} {output}"
rule spike_in_mpileup:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup"
    params:
        ref = '../reference.fasta/WT.ref.upper.fa'
    shell:
        "{SAMTOOLS} mpileup {input[0]} --reference {params.ref} --max-depth 10000000 -q 20 -Q 20 > {output} " 
rule parse_mpileup:
    input: 
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat"
    shell:
        "{PYTHON} ./program/parse-mpileup-V04.py -i {input} -o {output} -n 0"
rule bmat_plot:
    input: 
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat"
    output:
        "../all_plot/cutoff_{cutoff}.ext50/TargetSeq-{lib}_{sample}_cutoff_{cutoff}_indel.ext50.pdf"
    params:
        sgRNA_seq = 'TGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCG'
    shell:
        """
        if [[ `cat {input} |wc -l` -eq 1 ]]; then 
        echo "bmat is empty" 
        echo "will touch a empty file" 
        touch {output}
        touch {output}.empty.log
        echo "decrease the cutoff and try again" > {output}.empty.log
        else
        echo "bmat is ok!"
        echo "start to plot"
        #echo `cat {params.sgRNA_seq}`
        # sgRNA=`cat {params.sgRNA_seq}`
        sgRNA={params.sgRNA_seq}
        {PYTHON} ./program/plot-targetseq-bmat-V04.py -i {input} -o {output} --region_extend_length 50 --sgRNA $sgRNA
        fi
        """
