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

SAMPLES = ['HEK4-ABEABE8EACBECBEOUT-Detect1',
 'HEK4-ABEABE8EACBECBEOUT-Digenome34',
 'HEK4-ABEABE8EACBECBEOUT-N2',
 'HEK4-ABEABE8EACBEcbeOUT-C1',
 'HEK4-ABEABE8EACBEcbeout-B2',
 'HEK4-ABEABE8Eacbecbeout-D1',
 'HEK4-ABEABE8Eacbecbeout-D3',
 'HEK4-ABEABE8Eacbecbeout-D4',
 'HEK4-abeABE8EACBECBEout-A1',
 'HEK4-abeABE8EACBECBEout-A2',
 'HEK4-abeABE8EACBECBEout-A6',
 'HEK4-abeabe8eACBECBEout-A5',
 'HEK4-abeabe8eacbeCBEout-A8']
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
        expand('../reference.fasta/{sample}.bed',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.fa',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.tmp.fa',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.upper.fa',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.upper.fa.fai',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.upper.fa.amb',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.upper.fa.ann',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.upper.fa.bwt',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.upper.fa.pac',sample=SAMPLES),
        expand('../reference.fasta/{sample}.ref.upper.fa.sa',sample=SAMPLES),
        expand('../reference.fasta/{sample}.sgRNA.bed',sample=SAMPLES),
        expand('../reference.fasta/{sample}.sgRNA.fa',sample=SAMPLES),
        expand('../reference.fasta/{sample}.sgRNA.tmp.fa',sample=SAMPLES),
        expand('../reference.fasta/{sample}.sgRNA.upper.fa',sample=SAMPLES),
        expand('../reference.fasta/{sample}.sgRNA.upper.fa.seq',sample=SAMPLES)

rule getbed:
    output:
        temp('../reference.fasta/{sample}.bed'),
        temp('../reference.fasta/{sample}.sgRNA.bed')
    run:
        import os
        for tab in ['ABE']:
            path = './primer_table/{index}.txt'.format(index=tab)
            f = open(path,'r')
            ls = [x.split("\t") for x in f.readlines()][1:]
            # ls
            lib = "../reference.fasta"
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
        '../reference.fasta/{sample}.bed',
        '../reference.fasta/{sample}.sgRNA.bed'
    output:
        temp('../reference.fasta/{sample}.ref.tmp.fa'),
        temp('../reference.fasta/{sample}.sgRNA.tmp.fa')
    shell:
        """
        {BEDTOOLS} getfasta -nameOnly -s -bed {input[0]} -fi {GENOME_HG38} > {output[0]}
        {BEDTOOLS} getfasta -nameOnly -s -bed {input[1]} -fi {GENOME_HG38} > {output[1]}
        """
rule rm_charater_fa:
    input:
        '../reference.fasta/{sample}.ref.tmp.fa',
        '../reference.fasta/{sample}.sgRNA.tmp.fa'
    output:
        temp('../reference.fasta/{sample}.ref.fa'),
        temp('../reference.fasta/{sample}.sgRNA.fa')
    params:
        awk = """'FNR==1{print substr($1, 1, length($1)-3);} FNR==2{print;}'"""
    shell:
        """
        awk {params.awk} {input[0]} > {output[0]}
        awk {params.awk} {input[1]} > {output[1]}
        """
rule letter2LETTER:
    input:
        '../reference.fasta/{sample}.ref.fa',
        '../reference.fasta/{sample}.sgRNA.fa'
    output:
        '../reference.fasta/{sample}.ref.upper.fa',
        temp('../reference.fasta/{sample}.sgRNA.upper.fa')
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
        '../reference.fasta/{sample}.sgRNA.upper.fa'
    output:
        '../reference.fasta/{sample}.sgRNA.upper.fa.seq'
    shell:
        """
        sed '1d' {input} > {output}
        """
rule faidx:
    input:
        '../reference.fasta/{sample}.ref.upper.fa'
    output:
        '../reference.fasta/{sample}.ref.upper.fa.fai'
    shell:
        """
        {SAMTOOLS} faidx {input}
        """

rule bwa_idex:
    input:
        '../reference.fasta/{sample}.ref.upper.fa',
        '../reference.fasta/{sample}.ref.upper.fa.fai',
    output:
        '../reference.fasta/{sample}.ref.upper.fa.amb',
        '../reference.fasta/{sample}.ref.upper.fa.ann',
        '../reference.fasta/{sample}.ref.upper.fa.bwt',
        '../reference.fasta/{sample}.ref.upper.fa.pac',
        '../reference.fasta/{sample}.ref.upper.fa.sa'
    shell:
        """
        {BWA} index {input[0]}
        """