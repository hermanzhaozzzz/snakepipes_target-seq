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
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bedtools" # ok
SAMCLIP = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/samclip" # ok
PYTHON = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/python" # ok
GENOME_HG38 = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"



CUTOFF = ["3"]

LIBS = [
    '11',
    '12',
    '13'
]


SAMPLES = ['MTND6P4-N6',
 'VEGFA-OffNo-GDNo-08',
 'VEGFA-Off-Target-2',
 'VEGFA-Off-Target-3',
 'VEGFA-Off-Target-6',
 'VEGFA-Off-Target-7',
 'VEGFA-On-Target']

READ_IDX = ["1","2"]


rule all:
    input:
        expand("../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R{read_idx}_cutadapt.fq.gz",lib=LIBS,sample=SAMPLES,read_idx=READ_IDX,cutoff=CUTOFF),
        expand("../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF)
        


rule check_file:
    output:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq",
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R2.fastq"
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
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq",
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R2.fastq"
    output:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz"
    log:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_cutadapt.log"
    shell:
        """
        {CUTADAPT} -j 24 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
        """  

rule bwa_mapping:
    input:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz",
    output:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam"
    log:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.log"
    shell:
        """
        {BWA} mem {GENOME_HG38} {input[0]} {input[1]} -t 24 -M > {output} 2>{log}
        """
rule sam_to_bam:
    input:
        sam = "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam"
    output:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    shell:
        """
        {SAMTOOLS} view -h -f 1 -F 268 {input.sam} \
        | {SAMCLIP} --ref {GENOME_HG38} --max 3 --progress 0 \
        | awk 'function abs(v) {{return v < 0 ? -v : v}} $1~"@" || ($7 == "=" && abs($9) <= 600 ) {{print $0}}' \
        | {SAMTOOLS} view -hb > {output}
        """


rule samtools_sort_by_position:
    input:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    output:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    shell:
        "{SAMTOOLS} sort -O BAM -o {output} -T {output}.temp -@ 6 -m 2G {input}"


rule samtools_index:
    input:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    output:
        "../igv/TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai"
    shell:
        "{SAMTOOLS} index -@ 6 {input} {output}"