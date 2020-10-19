GENOME_HG38 = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"


BWA = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bwa"
SAMTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_cutadapt-STARmapping-FPKM-sortBAM/bin/samtools"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bedtools"


# 记得修改getbed中的txt地址
# 决定了txt table的数量！！！
# 不改会缺文件

SAMPLES = ['ABESite7-On',
 'EMX1Site2-On',
 'HEK2-On-Target',
 'HEK4-Guideseq-1',
 'HEK4-Guideseq-2',
 'HEK4-Guideseq-3',
 'HEK4-Guideseq-4',
 'HEK4-Guideseq-5',
 'HEK4-On-Target',
 'HEK5-On-Target',
 'MTND4P12-N4',
 'MTND5P11-N5.1',
 'MTND6P4-N6',
 'PP2-On-Target',
 'PP6-On-Target',
 'RNF2-ON-Target',
 'VEGFA-Off-Target-2',
 'VEGFA-Off-Target-3',
 'VEGFA-Off-Target-6',
 'VEGFA-Off-Target-7',
 'VEGFA-OffNo-GDNo-08',
 'VEGFA-OffNo-GDNo-12',
 'VEGFA-On-Target']

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
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.sa',sample=SAMPLES)
rule getbed:
    output:
        '../TargetSeq_BED_sample_lib/{sample}.bed'
    run:
#         import os
# #         path = "./primer_table/primer_info.txt"
#         path = 
#         f = open(path,'r')
#         ls = [x.split("\t") for x in f.readlines()][1:]
#         # ls
#         lib = "../TargetSeq_BED_sample_lib"
#         try:
#             os.makedirs(lib)
#         except:
#             pass
#         for ls_line in ls:
#         #     print(ls_line)
#             flag = ls_line[3]
#             chr_ = ls_line[0]
#             tss = str(int(ls_line[1])-1)
#             tes = ls_line[2]
#             strand = ls_line[8]
#             with open(os.path.join(lib,flag)+".bed",'w') as bed:
#                 bed.write('{chr}\t{tss}\t{tes}\t{flag}\t0\t{strand}\n'.format(chr=chr_,tss=tss,tes=tes,flag=flag,strand=strand))
        import os
        # 1.txt 2.txt 3.txt.......7.txt
        for tab in ['5-6-7','8-9-10','11-12-13','14-15','duigou-jiantou','quanquan-sanjiao']:
            path = './primer_table/{index}.txt'.format(index=tab)
            f = open(path,'r')
            ls = [x.split("\t") for x in f.readlines()][1:]
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
                tss = str(int(ls_line[1])-1)
                tes = ls_line[2]
                strand = ls_line[8]
                with open(os.path.join(lib,flag)+".bed",'w') as bed:
                    bed.write('{chr}\t{tss}\t{tes}\t{flag}\t0\t{strand}\n'.format(chr=chr_.strip(),tss=tss.strip(),tes=tes.strip(),flag=flag.strip(),strand=strand.strip()))
rule getfasta:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.bed'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.tmp.fa'
    shell:
        """
        {BEDTOOLS} getfasta -nameOnly -s -bed {input} -fi {GENOME_HG38} > {output}
        """
rule rm_charater_fa:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.tmp.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa'
    params:
        awk = """'FNR==1{print substr($1, 1, length($1)-3);} FNR==2{print;}'"""
    shell:
        """
        awk {params.awk} {input} > {output}
        """
rule letter2LETTER:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
    params:
        awk = """'{print toupper($0)}'"""
    shell:
        """
        cat {input} | awk {params.awk} > {output}
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
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.amb',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.ann',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.bwt',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.pac',
        '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa.sa'        
    shell:
        """
        {BWA} index {input}
        """