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
    'ND1-on-target_ND1-L1397N-rep1-0301',
    'ND1-on-target_ND1-L1397N-rep2-0301',
    'ND1-on-target_untreated-rep1-0301',
    'ND1-on-target_untreated-rep2-0301',
    'ND4-on-target_ND4-2000-12-rep1-0301',
    'ND4-on-target_ND4-2000-12-rep2-0301',
    'ND4-on-target_ND4-L1333N-rep1-0301',
    'ND4-on-target_ND4-L1333N-rep2-0301',
    'ND4-on-target_ND4-L1397C-rep1-0301',
    'ND4-on-target_ND4-L1397C-rep2-0301',
    'ND4-on-target_untreated-rep1-0301',
    'ND4-on-target_untreated-rep2-0301',
    'ND5.1-on-target_DddA-NES-ND5.1-rep1-0301',
    'ND5.1-on-target_DddA-NES-ND5.1-rep2-0301',
    'ND5.1-on-target_ND5.1-1333-DddIA-rep1-0301',
    'ND5.1-on-target_ND5.1-1333-DddIA-rep2-0301',
    'ND5.1-on-target_ND5.1-1333-mapKNES-rep1-0301',
    'ND5.1-on-target_ND5.1-1333-rep1-0301',
    'ND5.1-on-target_ND5.1-1333-rep2-0301',
    'ND5.1-on-target_ND5.1-2000-12-rep1-0301',
    'ND5.1-on-target_ND5.1-2000-12-rep2-0301',
    'ND5.1-on-target_TALE-NES-ND5.1-rep1-0301',
    'ND5.1-on-target_TALE-NES-ND5.1-rep2-0301',
    'ND5.1-on-target_UGI-NES-ND5.1-rep1-0301',
    'ND5.1-on-target_UGI-NES-ND5.1-rep2-0301',
    'ND5.1-on-target_WT-ND5.1-rep1-0301',
    'ND5.1-on-target_WT-ND5.1-rep2-0301',
    'ND5.1-on-target_untreated-rep1-0301',
    'ND5.1-on-target_untreated-rep2-0301',
    'ND5.3-on-target_ND5.3-L1397C-rep1-0301',
    'ND5.3-on-target_ND5.3-L1397C-rep2-0301',
    'ND5.3-on-target_untreated-rep1-0301',
    'ND5.3-on-target_untreated-rep2-0301',
    'ND6-on-target_DddA-NES-ND6-rep1',
    'ND6-on-target_DddA-NES-ND6-rep2',
    'ND6-on-target_DddIA-1-0.25-ND6-rep1',
    'ND6-on-target_DddIA-1-0.25-ND6-rep2',
    'ND6-on-target_DddIA-1-0.5-ND6-rep1',
    'ND6-on-target_DddIA-1-0.5-ND6-rep2',
    'ND6-on-target_DddIA-1-1-ND6-rep1',
    'ND6-on-target_DddIA-1-1-ND6-rep2',
    'ND6-on-target_DddIA-1-1.2-ND6-rep1',
    'ND6-on-target_DddIA-1-1.2-ND6-rep2',
    'ND6-on-target_DddIA-1-1.5-ND6-rep1',
    'ND6-on-target_DddIA-1-1.5-ND6-rep2',
    'ND6-on-target_G1309A-ND6-rep1',
    'ND6-on-target_G1309A-ND6-rep2',
    'ND6-on-target_N1308A-ND6-rep1',
    'ND6-on-target_N1308A-ND6-rep2',
    'ND6-on-target_N1367A-ND6-rep1',
    'ND6-on-target_N1367A-ND6-rep2',
    'ND6-on-target_N1368A-ND6-rep1',
    'ND6-on-target_N1368A-ND6-rep2-0301',
    'ND6-on-target_N6-2000-12-rep1-0301',
    'ND6-on-target_N6-2000-12-rep2-0301',
    'ND6-on-target_N6-LTX-12-rep1-0301',
    'ND6-on-target_N6-LTX-12-rep2-0301',
    'ND6-on-target_ND6-1397-DddIA-rep1-0301',
    'ND6-on-target_ND6-1397-DddIA-rep2-0301',
    'ND6-on-target_ND6-1397-hivNES-Q1310A-rep1-0301',
    'ND6-on-target_ND6-1397-hivNES-Q1310A-rep2-0301',
    'ND6-on-target_ND6-1397-hivNES-rep1-0301',
    'ND6-on-target_ND6-1397-hivNES-rep2-0301',
    'ND6-on-target_ND6-1397-rep1-0301',
    'ND6-on-target_ND6-1397-rep2-0301',
    'ND6-on-target_Q1310A-ND6-rep1',
    'ND6-on-target_Q1310A-ND6-rep2',
    'ND6-on-target_TALE-NES-ND6-rep1',
    'ND6-on-target_TALE-NES-ND6-rep2',
    'ND6-on-target_UGI-NES-ND6-rep1',
    'ND6-on-target_UGI-NES-ND6-rep2',
    'ND6-on-target_WT-ND6-rep1',
    'ND6-on-target_WT-ND6-rep2',
    'ND6-on-target_untreated-rep1',
    'ND6-on-target_untreated-rep2',
    # 'test'
]

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
PYTHON = "/lustre1/chengqiyi_pkuhpc/zhaohn/miniconda3/envs/snakepipes_py27/bin/python"
# abyss
# PYTHON = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq/bin/python"

dt = {
    "ND1-on-target": "ACTCAATCCTCTGATC",
    "ND4-on-target": "CCTGATCAAATATC",
    "ND5.1-on-target": "ACCTTTCCTCACAGGT",
    "ND5.3-on-target": "CTACTCATCTTCCTAATT",
    "ND6-on-target": "CCTCAGGATACTCCT"
}

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
        expand("../../1.5generation/mapping/{sample}_R1_cutadapt.fq.gz",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.bam",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.bam.bai",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.mpileup",sample=SAMPLES),
        expand("../../1.5generation/mapping/{sample}_bwa_sort.bmat",sample=SAMPLES),
        expand("../../all_plot/1.5generation/{sample}.ext50.pdf",sample=SAMPLES),
rule cutadapt:
    input:
        "../../fastq/{sample}_R1.fastq"
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
        ref = lambda wildcards, input: f"../../reference.fasta/{input[0].split('mapping/')[1].split('_')[0]}.ref.upper.fa",
    log:
        "../../1.5generation/mapping/{sample}_bwa.log"
    shell:
        """
        {BWA} mem {params.ref} {input} -t {THREADS} -M > {output} 2>{log}
        """
rule sam_to_bam:
    input:
        sam = "../../1.5generation/mapping/{sample}_bwa.sam",
    output:
        temp("../../1.5generation/mapping/{sample}_bwa.bam")
    shell:
        """
        {SAMTOOLS} view -hb {input.sam}> {output}
        """
rule samtools_sort_by_position:
    input:
        "../../1.5generation/mapping/{sample}_bwa.bam"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam"
    shell:
        "{SAMTOOLS} sort -O BAM -o {output} -T {output}.temp -@ {THREADS} -m 2G {input}"
rule samtools_index:
    input:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam.bai"
    shell:
        "{SAMTOOLS} index -@ {THREADS} {input} {output}"
rule spike_in_mpileup:
    input:
        "../../1.5generation/mapping/{sample}_bwa_sort.bam",
        "../../1.5generation/mapping/{sample}_bwa_sort.bam.bai"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.mpileup"
    params:
        ref = lambda wildcards, input: f"../../reference.fasta/{input[0].split('mapping/')[1].split('_')[0]}.ref.upper.fa",
    shell:
        "{SAMTOOLS} mpileup {input[0]} --reference {params.ref} --max-depth 10000000 -q 20 -Q 20 > {output} " 
rule parse_mpileup:
    input: 
        "../../1.5generation/mapping/{sample}_bwa_sort.mpileup"
    output:
        "../../1.5generation/mapping/{sample}_bwa_sort.bmat"
    shell:
        "{PYTHON} ../program/parse-mpileup-V04.py -i {input} -o {output} -n 0"
rule bmat_plot:
    input: 
        "../../1.5generation/mapping/{sample}_bwa_sort.bmat"
    output:
        "../../all_plot/1.5generation/{sample}.ext50.pdf"
    params:
        sgRNA_seq = lambda wildcards, input: dt[input[0].split("mapping/")[1].split("_")[0]],
    shell:
        """
        if [[ `cat {input} |wc -l` -eq 1 ]]; then 
        echo "bmat is empty" 
        echo "will touch a empty file" 
        touch {output}
        touch {output}.empty.log
        else
        echo "bmat is ok!"
        echo "start to plot"
        #echo `cat {params.sgRNA_seq}`
        # sgRNA=`cat {params.sgRNA_seq}`
        sgRNA={params.sgRNA_seq}
        {PYTHON} ../program/plot-targetseq-bmat-V04.py -i {input} -o {output} --region_extend_length 50 --sgRNA $sgRNA
        fi
        """
