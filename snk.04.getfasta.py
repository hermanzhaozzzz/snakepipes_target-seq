GENOME_HG38 = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"


BWA = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bwa"
SAMTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_cutadapt-STARmapping-FPKM-sortBAM/bin/samtools"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bedtools"


# 记得修改getbed中的txt地址
# 决定了txt table的数量！！！
# 不改会缺文件

SAMPLES = ['ND516-share-4',
 'ND516-share-5',
 'ND516-share-9',
 'ND6-only-11',
 'ND6-only-2',
 'ND6-only-6']

rule all:
    input:
        expand('../TargetSeq_BED_sample_lib/{sample}.bed',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.tmp.fa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.fai',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.amb',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.ann',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.bwt',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.pac',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.sa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.sgRNA.bed',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.sgRNA.fa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.sgRNA.tmp.fa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.sgRNA.upper.fa',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.sgRNA.upper.fa.seq',sample=SAMPLES)

rule getbed:
    output:
        '../TargetSeq_BED_sample_lib/{sample}.bed',
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.bed'
    run:
        import os
        # 1.txt 2.txt 3.txt.......7.txt
        for tab in ['sub1']:
            path = './primer_table/{index}.txt'.format(index=tab)
            f = open(path,'r')
            ls = [x.split("\t") for x in f.readlines()][1:]
            # ls
            lib = "../TargetSeq_BED_sample_lib"
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
        '../TargetSeq_BED_sample_lib/{sample}.bed',
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.bed'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.tmp.fa',
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.tmp.fa'
    shell:
        """
        {BEDTOOLS} getfasta -nameOnly -s -bed {input[0]} -fi {GENOME_HG38} > {output[0]}
        {BEDTOOLS} getfasta -nameOnly -s -bed {input[1]} -fi {GENOME_HG38} > {output[1]}
        """
rule rm_charater_fa:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.tmp.fa',
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.tmp.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa',
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.fa'
    params:
        awk = """'FNR==1{print substr($1, 1, length($1)-3);} FNR==2{print;}'"""
    shell:
        """
        awk {params.awk} {input[0]} > {output[0]}
        awk {params.awk} {input[1]} > {output[1]}
        """
rule letter2LETTER:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa',
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa',
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.upper.fa'
    params:
        awk = """'{print toupper($0)}'"""
    shell:
        """
        cat {input[0]} | awk {params.awk} > {output[0]}
        cat {input[1]} | awk {params.awk} > {output[1]}
        """
rule get_sgRNAseq:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.upper.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.sgRNA.upper.fa.seq'
    shell:
        """
        sed '1d' {input} > {output}
        """
rule faidx:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.fai'
    shell:
        """
        {SAMTOOLS} faidx {input}
        """

rule bwa_idex:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.fai',
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.amb',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.ann',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.bwt',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.pac',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.sa'
    shell:
        """
        {BWA} index {input[0]}
        """