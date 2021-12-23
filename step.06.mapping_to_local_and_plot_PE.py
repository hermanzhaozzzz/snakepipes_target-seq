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


CUTOFF = [
    "3",
    # "5","10"
]

LIBS = ['untreated-rep1', 'untreated-rep2']


SAMPLES = ['share-new-1',
 'share-new-2',
 'share-new-3',
 'share-new-4',
 'share-new-5',
 'share-new-6',
 'share-new-7',
 'share-new-8',
 'share-new-9',
 'share-new-10',
 'share-new-11',
 'share-new-12',
 'share-new-13',
 'share-new-14',
 'share-new-15',
 'share-new-16',
 'share-new-17',
 'share-new-18',
 # 'share-new-19', # extend
 'share-new-20',
 'share-new-21',
 'share-new-22',
 'share-new-23',
 'share-new-24',
 'share-new-25',
 'share-new-26',
 'share-new-27',
 'share-new-28',
 'share-new-29',
 'share-new-30',
 'share-new-31',
 'share-new-32',
 'ND4-new-only1',
 'ND4-new-only2',
 'ND4-new-only3',
 'ND4-new-only4',
 'ND4-new-only5',
 'ND4-new-only6',
 'ND4-new-only7',
 'ND4-new-only8',
 'ND4-new-only9',
 'ND4-new-only10',
 'ND4-new-only11',
 'ND4-new-only12',
 'ND4-new-only13',
 'ND5.1-new-only1',
 'ND5.1-new-only2',
 'ND5.1-new-only3',
 'ND5.1-new-only4',
 'ND5.1-new-only5',
 'ND5.1-new-only6',
 'ND5.1-new-only7',
 'ND5.1-new-only8',
 'ND5.1-new-only9',
 'ND5.1-new-only10',
 'ND5.1-new-only11',
 'ND5.1-new-only12',
 'ND5.1-new-only13',
 'ND5.1-new-only14',
 'ND5.1-new-only15',
 'ND5.1-new-only16',
 'ND5.1-new-only17',
 'ND5.1-new-only18',
 'ND5.1-new-only19',
 'ND5.1-new-only20',
 'ND5.1-new-only21',
 'ND5.1-new-only22',
 'ND5.1-new-only23',
 'ND5.1-new-only24',
 'ND5.1-new-only25',
 'ND6-new-only1',
 'ND6-new-only2',
 'ND6-new-only3',
 'ND6-new-only4',
 'ND6-new-only5',
 'ND6-new-only6',
 'ND6-new-only7',
 'ND6-new-only8',
 'ND6-new-only9',
 'ND6-new-only10',
 'ND6-new-only11',
 'ND6-new-only12',
 'ND6-new-only13',
 'ND6-new-only14',
 'ND6-new-only15',
 'ND6-new-only16',
 'ND6-new-only17']

READ_IDX = ["1","2"]

# defult_sgRNA_dict_for_plot = {
#     "EMX1": "GAGTCCGAGCAGAAGAAGAAGGG",
# }

# defult_sgRNA_dict_for_plot = {
#     "ABESite7":"GAATACTAAGCATAGACTCC", # 这里指的ABE-Site-7-On-Target
#     "EMX1": "GAGTCCGAGCAGAAGAAGAAGGG",
#     'EMX1Site2': "GTATTCACCTGAAAGTGTGC", 
#     'HEK2':'GAACACAAAGCATAGACTGC',
#     'HEK3': "GGCCCAGACTGAGCACGTGATGG", 
#     "HEK4":"GGCACTGCGGCTGGAGGTGGGGG", 
#     "HEK5":'CTGGCCTGGGTCAATCCTTG',
#     "MTND4P12":"TGCTAGTAACCACATTCTCCTGATCAAATATCACTCTCCTACTTACAGGA",
#     "MTND5P11":"TAGCATTGGCAGGAATACCCTTCCTCACAGGTTTCTACTCCAAAGA",
#     "MTND6P4":"TGACCCCCATGCCTCAGGATACTCCTCAATAGCCACCG",
#     "PP2":"GGCACTCGGGGGCGAGAGGA",
#     "PP6":"GGGGCTCAACATCGGAAGAG",
#     "RNF2":"GTCATCTTAGTCATTACCTGAGG",
#     'VEGFA': "GACCCCCTCCACCCCGCCTCCGG",
#     "FANCF": "GGAATCCCTTCTGCAGCACC",
#     "ND6":"TGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCG",
#     "ND5.1":"TAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGA",
#     "ND4":"TGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGA",
#     "ND516":"TAGCATTAGCAGGAATACCTTTC7CTCACAGGTTTCTACTCCAAAGA", # run ND5.1 genome use ND5.1 on-target sequence
#     "ND516":"TGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCG", # run ND6 genome use ND6 on-target sequence
# }


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
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R{read_idx}_cutadapt.fq.gz",lib=LIBS,sample=SAMPLES,read_idx=READ_IDX,cutoff=CUTOFF),
        expand('../reference.fasta/{sample}.ref.upper.fa',sample=SAMPLES),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../all_plot/cutoff_{cutoff}.ext50/TargetSeq-{lib}_{sample}_cutoff_{cutoff}_indel.ext50.pdf",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand('../reference.fasta/{sample}.sgRNA.upper.fa.seq',sample=SAMPLES)
rule check_file:
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq",
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R2.fastq"
    run:
        try:
            open(output[0],'r')
        except IOError:
            open(output[0],'w')
        try:
            open(output[1],'r')
        except IOError:
            open(output[1],'w')
            
rule cutadapt:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq",
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R2.fastq"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz"
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_cutadapt.log"
    shell:
        """
        {CUTADAPT} -j {THREADS} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
        """
rule bwa_mapping:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz",
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam"
    params:
        ref = '../reference.fasta/{sample}.ref.upper.fa'
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.log"
    shell:
        """
        {BWA} mem {params.ref} {input[0]} {input[1]} -t {THREADS} -M > {output} 2>{log}
        """
# rule bowtie_build_index:
#     input:
#         '../reference.fasta/{sample}.ref.upper.fa'
#     output:
#         '../reference.fasta/{sample}.ref.upper.fa.bowtie_index'# 形式输出1个即可
#     log:
#         '../reference.fasta/{sample}.ref.upper.fa.bowtie-build.log'
#     shell:
#         """
#         touch {output}
#         echo "a formal output test" > {output}
#         {BOWTIE_BUILD} {input} {output} 2>{log}
#         """
# rule bowtie1_mapping:
#     input:
#         "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
#         "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz",
#         '../reference.fasta/{sample}.ref.upper.fa.bowtie_index'
#     output:
#         "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam"
#     log:
#         "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.log"
#     shell:
#         """
#         {BOWTIE} -x {input[2]} -1 {input[0]} -2 {input[1]} -p {THREADS} -S {output} 2>{log}
#         """

##
rule sam_to_bam:
    input:
        sam = "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam",
        ref = '../reference.fasta/{sample}.ref.upper.fa'
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    shell:
        """
        {SAMTOOLS} view -h -f 1 -F 268 {input.sam} \
        | {SAMCLIP} --ref {input.ref} --max 3 --progress 0 \
        | awk 'function abs(v) {{return v < 0 ? -v : v}} $1~"@" || ($7 == "=" && abs($9) <= 600 ) {{print $0}}' \
        | {SAMTOOLS} view -hb > {output}
        """


rule samtools_sort_by_position:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    shell:
        "{SAMTOOLS} sort -O BAM -o {output} -T {output}.temp -@ {THREADS} {input}"


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
        ref = '../reference.fasta/{sample}.ref.upper.fa'
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
#         sgRNA_seq = lambda wildcards, input: defult_sgRNA_dict_for_plot[input[0].split("/")[4].split("-")[0]],
        sgRNA_seq = '../reference.fasta/{sample}.sgRNA.upper.fa.seq'
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
        sgRNA=`cat {params.sgRNA_seq}`
        # sgRNA={params.sgRNA_seq}
        {PYTHON} ./program/plot-targetseq-bmat-V04.py -i {input} -o {output} --region_extend_length 50 --sgRNA $sgRNA
        fi
        """
