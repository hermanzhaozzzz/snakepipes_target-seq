GENOME_HG38 = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
BWA = "/home/zhaohuanan/miniconda3/bin/bwa"
SAMTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_cutadapt-STARmapping-FPKM-sortBAM/bin/samtools"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/py27/bin/bedtools"


SAMPLES = [
    "HK4-mcf-1",
    "HK4-mcf-2",
    "HK4-mcf-3",
    "HK4-mcf-4",
    "HK4-mcf-5",
    "HK4-mcf-6",
    "HK4-mcf-7",
    "HK4-mcf-8",
    "HK4-mcf-9",
    "HK4-mcf-10",
    "HK4-mcf-11",
    "HK4-mcf-12",
    "HK4-mcf-13",
    "HK4-mcf-14",
    "HK4-mcf-15"
]

rule all:
    input:
        expand('../TargetSeq_BED_sample_lib/{sample}.bed',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa.fai',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa.amb',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa.ann',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa.bwt',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa.pac',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa.sa',sample=SAMPLES)
rule getbed:
    output:
        '../TargetSeq_BED_sample_lib/{sample}.bed'
    run:
        import os
        path = "./primer_table/primer_info.txt"
        f = open(path,'r')
        ls = [x.split("\t")[:4] for x in f.readlines()][1:]
        # ls
        lib = "../TargetSeq_BED_sample_lib"
        try:
            os.makedirs(lib)
        except:
            pass
        for ls_line in ls:
        #     print(ls_line)
            flag = ls_line[3]
            chr_ = ls_line[0]
            tss = ls_line[1]
            tes = ls_line[2]
            with open(os.path.join(lib,flag)+".bed",'w') as bed:
                bed.write('{chr}\t{tss}\t{tes}\n'.format(chr=chr_,tss=tss,tes=tes))
rule getfasta:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.bed'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa'
    shell:
        """
        {BEDTOOLS} getfasta -bed {input} -fi {GENOME_HG38} > {output}
        """

rule faidx:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa.fai'
    shell:
        """
        {SAMTOOLS} faidx {input}
        """
rule bwa_idex:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa.amb',
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa.ann',
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa.bwt',
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa.pac',
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa.sa'        
    shell:
        """
        {BWA} index {input}
        """