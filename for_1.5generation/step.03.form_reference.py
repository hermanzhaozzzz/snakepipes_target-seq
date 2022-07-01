import os

# ------------------------------------------------------------------->>>>>>>>>>
# Author: Herman ZHAO
# Email: hermanzhaozzzz@gmail.com
# Update
#     date: 2021-08-03 log: format annotation and software path
# ------------------------------------------------------------------->>>>>>>>>>

# ------------------------------------------------------------------->>>>>>>>>>
# region_id INFO
# ------------------------------------------------------------------->>>>>>>>>>

REGION_ID = [
    # "ND1-on-target",
    # "ND4-on-target",
    # "ND5.1-on-target",
    # "ND5.3-on-target",
    "ND6-on-target",
]
# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
THREADS = '20'

# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
PRIMER_INFO = "primer.txt"
GENOME_HG38 = "/lustre1/chengqiyi_pkuhpc/zhaohn/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# polaris
PYTHON = "/lustre1/chengqiyi_pkuhpc/zhaohn/miniconda3/envs/snakepipes_py27/bin/python"
# abyss
# PYTHON = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq/bin/python"

# get the application path
with os.popen("which bwa") as path:
    BWA = path.read().strip()
    print('PATH bwa:', BWA)
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
    print('PATH samtools:', SAMTOOLS)
with os.popen("which bedtools") as path:
    BEDTOOLS = path.read().strip()
    print('PATH bedtools:', BEDTOOLS)
if "" in [BWA, SAMTOOLS, BEDTOOLS]:
    raise ValueError("The necessary software path is missing!!!")
else:
    print("pass software check...\n")

rule all:
    input:
        expand('../../reference.fasta/{region_id}.bed',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.fa',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.tmp.fa',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.upper.fa',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.upper.fa.fai',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.upper.fa.amb',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.upper.fa.ann',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.upper.fa.bwt',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.upper.fa.pac',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.ref.upper.fa.sa',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.sgRNA.bed',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.sgRNA.fa',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.sgRNA.tmp.fa',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.sgRNA.upper.fa',region_id=REGION_ID),
        expand('../../reference.fasta/{region_id}.sgRNA.upper.fa.seq',region_id=REGION_ID)

rule getbed:
    output:
        temp('../../reference.fasta/{region_id}.bed'),
        temp('../../reference.fasta/{region_id}.sgRNA.bed')
    run:
        import os
        for tab in ['primer']:
            path = '{index}.txt'.format(index=tab)
            f = open(path,'r')
            ls = [x.split("\t") for x in f.readlines()][1:]
            # ls
            lib = "../../reference.fasta"
            try:
                os.makedirs(lib)
            except:
                pass
            # form ref bed
            # form sgRNA bed 2
            for ls_line in ls:
                print(ls_line)
                flag = ls_line[3]
                chr_ = ls_line[0]
                tss = str(int(ls_line[1])-1)
                tes = ls_line[2]
                strand = ls_line[8]
                print('{chr}\t{tss}\t{tes}\t{flag}\t0\t{strand}\n'.format(chr=chr_.strip(),tss=tss.strip(),tes=tes.strip(),flag=flag.strip(),strand=strand.strip()))
                flag2 = ls_line[3]
                chr_2 = ls_line[10]
                tss2 = str(int(ls_line[11])-1)
                tes2 = ls_line[12]
                if tss2 < tss:
                    tss2 = tss
                if tes2 > tes:
                    tes2 = tes
                strand2 = ls_line[8]
                print('{chr}\t{tss}\t{tes}\t{flag}\t0\t{strand}\n'.format(chr=chr_2.strip(),tss=tss2.strip(),tes=tes2.strip(),flag=flag2.strip(),strand=strand2.strip()))
                with open(os.path.join(lib,flag)+".bed",'w') as bed:
                    bed.write('{chr}\t{tss}\t{tes}\t{flag}\t0\t{strand}\n'.format(chr=chr_.strip(),tss=tss.strip(),tes=tes.strip(),flag=flag.strip(),strand=strand.strip()))
                with open(os.path.join(lib,flag)+".sgRNA.bed",'w') as bed:
                    bed.write('{chr}\t{tss}\t{tes}\t{flag}\t0\t{strand}\n'.format(chr=chr_2.strip(),tss=tss2.strip(),tes=tes2.strip(),flag=flag2.strip(),strand=strand2.strip()))
rule getfasta:
    input:
        '../../reference.fasta/{region_id}.bed',
        '../../reference.fasta/{region_id}.sgRNA.bed'
    output:
        temp('../../reference.fasta/{region_id}.ref.tmp.fa'),
        temp('../../reference.fasta/{region_id}.sgRNA.tmp.fa')
    shell:
        """
        {BEDTOOLS} getfasta -nameOnly -s -bed {input[0]} -fi {GENOME_HG38} > {output[0]}
        {BEDTOOLS} getfasta -nameOnly -s -bed {input[1]} -fi {GENOME_HG38} > {output[1]}
        """
rule rm_charater_fa:
    input:
        '../../reference.fasta/{region_id}.ref.tmp.fa',
        '../../reference.fasta/{region_id}.sgRNA.tmp.fa'
    output:
        temp('../../reference.fasta/{region_id}.ref.fa'),
        temp('../../reference.fasta/{region_id}.sgRNA.fa')
    params:
        awk = """'FNR==1{print substr($1, 1, length($1)-3);} FNR==2{print;}'"""
    shell:
        """
        awk {params.awk} {input[0]} > {output[0]}
        awk {params.awk} {input[1]} > {output[1]}
        """
rule letter2LETTER:
    input:
        '../../reference.fasta/{region_id}.ref.fa',
        '../../reference.fasta/{region_id}.sgRNA.fa'
    output:
        '../../reference.fasta/{region_id}.ref.upper.fa',
        temp('../../reference.fasta/{region_id}.sgRNA.upper.fa')
    params:
        awk = """'FNR==1{print $0;} FNR==2{print toupper($0);}'"""
#         awk = """'{print toupper($0)}'"""
    shell:
        """
        cat {input[0]} | awk {params.awk} > {output[0]}
        cat {input[1]} | awk {params.awk} > {output[1]}
        """
rule get_sgRNAseq:
    input:
        '../../reference.fasta/{region_id}.sgRNA.upper.fa'
    output:
        '../../reference.fasta/{region_id}.sgRNA.upper.fa.seq'
    shell:
        """
        sed '1d' {input} > {output}
        """
rule faidx:
    input:
        '../../reference.fasta/{region_id}.ref.upper.fa'
    output:
        '../../reference.fasta/{region_id}.ref.upper.fa.fai'
    shell:
        """
        {SAMTOOLS} faidx {input}
        """

rule bwa_idex:
    input:
        '../../reference.fasta/{region_id}.ref.upper.fa',
        '../../reference.fasta/{region_id}.ref.upper.fa.fai',
    output:
        '../../reference.fasta/{region_id}.ref.upper.fa.amb',
        '../../reference.fasta/{region_id}.ref.upper.fa.ann',
        '../../reference.fasta/{region_id}.ref.upper.fa.bwt',
        '../../reference.fasta/{region_id}.ref.upper.fa.pac',
        '../../reference.fasta/{region_id}.ref.upper.fa.sa'
    shell:
        """
        {BWA} index {input[0]}
        """