# _*_ coding: UTF-8 _*_

########################################################################
# MENG Howard
# 2020-07-07
# run demultiplex target seq [ APO Mut site ]
########################################################################
# run on abyss
# /home/menghaowei/menghw_HD/BE_project/20.target_seq_all/16.APOmut_screen.20200707
CUTADAPT = "/home/zhaohuanan/miniconda3/envs/cutadapt/bin/cutadapt"
BWA = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bwa"
SAMTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_cutadapt-STARmapping-FPKM-sortBAM/bin/samtools"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bedtools"
SAMCLIP = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/samclip"
PYTHON = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/python"





CUTOFF = ["3"]

LIBS = [
    'EH-1',
    'EH-2'
]

SAMPLES = ['EMX1-notOFF-01',
 'EMX1-notOFF-02',
 'EMX1-notOFF-03',
 'HEK4-Guideseq-3',
 'HEK4-Guideseq-4',
 'HEK4-Guideseq-6',
 'HEK4-Guideseq-3-rev',
 'HEK4-Guideseq-4-rev',
 'HEK4-Guideseq-6-rev']

READ_IDX = ["1","2"]

defult_sgRNA_dict_for_plot = {
    'VEGFA': "GACCCCCTCCACCCCGCCTCCGG", 
    'EMX1': "GAGTCCGAGCAGAAGAAGAAGGG", 
    'HEK3': "GGCCCAGACTGAGCACGTGATGG", 
    "HEK4":"GGCACTGCGGCTGGAGGTGGGGG", 
    "RNF2":"GTCATCTTAGTCATTACCTGAGG"
}
# geneSite = SAMPLES[0].split('-')[0]
# sgRNA_seq = defult_sgRNA_dict_for_plot[geneSite]
# print("The aim GENE site is:", geneSite, sgRNA_seq)

# print("Start snake rules")

rule all:
    input:
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R{read_idx}_cutadapt.fq.gz",lib=LIBS,sample=SAMPLES,read_idx=READ_IDX,cutoff=CUTOFF),

def abstract_gene_site(sample):
    geneSite = sample.split('-')[0]
    sgRNA_seq = defult_sgRNA_dict_for_plot[geneSite]
    return sgRNA_seq,geneSite
rule cutadapt:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq",
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R2.fastq"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz"
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_cutadapt.log"
    params:
        test = lambda wildcards, input: defult_sgRNA_dict_for_plot[input[0].split("/")[4].split("-")[0]]
    shell:
        #"""
        #srun -T 24 \
        """
        {CUTADAPT} -j 24 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
        echo this is {params.test}!!
        """  