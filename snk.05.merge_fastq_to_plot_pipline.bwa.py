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




CUTOFF = ["3","5","10"]

LIBS = ['1341P',
 'DddA-nes',
 'DddIA-1-1-0.5',
 'DddIA-1-1-1',
 'DddIA-1-1-1.5',
 'ND6-WT',
 'TALE-nes',
 'UGI-nes']

SAMPLES = ['ND516-share-4',
 'ND516-share-5',
 'ND516-share-9',
 'ND6-only-2',
 'ND6-only-6',
 'ND6-only-11']

READ_IDX = ["1","2"]

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
#     "ND5.1":"TAGCATTAGCAGGAATACCTTTC7CTCACAGGTTTCTACTCCAAAGA",
#     "ND4":"TGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGA",
# #     "ND516":"TAGCATTAGCAGGAATACCTTTC7CTCACAGGTTTCTACTCCAAAGA", # run ND5.1 genome use ND5.1 on-target sequence
#     "ND516":"TGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCG", # run ND6 genome use ND6 on-target sequence
# }

rule all:
    input:
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R{read_idx}_cutadapt.fq.gz",lib=LIBS,sample=SAMPLES,read_idx=READ_IDX,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../all_plot/cutoff_{cutoff}.ext50/TargetSeq-{lib}_{sample}_cutoff_{cutoff}_indel.ext50.pdf",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand('../TargetSeq_BED_sample_lib/{sample}.sgRNA.upper.fa.seq',sample=SAMPLES)
        

# wildcard_constraints:
#     pass


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
        #"""
        #srun -T 24 \
        """
        {CUTADAPT} -j 24 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
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
        ref = '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.log"
    shell:
        """
        {BWA} mem {params.ref} {input[0]} {input[1]} -t 24 -M > {output} 2>{log}
        """
rule sam_to_bam:
    input:
        sam = "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam",
        ref = '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
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
        "{SAMTOOLS} sort -O BAM -o {output} -T {output}.temp -@ 6 -m 2G {input}"


rule samtools_index:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai"
    shell:
        "{SAMTOOLS} index -@ 6 {input} {output}"


rule spike_in_mpileup:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup"
    params:
        ref = '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
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
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat",
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.upper.fa.seq'
    output:
        "../all_plot/cutoff_{cutoff}.ext50/TargetSeq-{lib}_{sample}_cutoff_{cutoff}_indel.ext50.pdf"
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
        echo `cat {input[1]}`
        sgRNA=`cat {input[1]}`
        {PYTHON} ./program/plot-targetseq-bmat-V04.py -i {input[0]} -o {output} --region_extend_length 50 --sgRNA $sgRNA
        fi
        """