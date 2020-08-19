GENOME_HG38 = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"


BWA = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bwa"
SAMTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_cutadapt-STARmapping-FPKM-sortBAM/bin/samtools"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bedtools"


# 记得修改getbed中的txt地址
# 决定了txt table的数量！！！
# 不改会缺文件

SAMPLES = ['EMX1-Dis-1',
 'EMX1-Dis-2',
 'EMX1-Dis-3',
 'EMX1-guide-1',
 'EMX1-guide-10',
 'EMX1-guide-13',
 'EMX1-guide-2',
 'EMX1-guide-4',
 'EMX1-notOFF-01',
 'EMX1-notOFF-02',
 'EMX1-notOFF-03',
 'EMX1-on-target',
 'HEK4-Dis-1',
 'HEK4-Dis-5',
 'HEK4-Guideseq-3',
 'HEK4-Guideseq-4',
 'HEK4-Guideseq-6',
 'VEGFA-Dis-1',
 'VEGFA-Dis-10',
 'VEGFA-Dis-11',
 'VEGFA-Dis-12',
 'VEGFA-Dis-13',
 'VEGFA-Dis-14',
 'VEGFA-Dis-2',
 'VEGFA-Dis-3',
 'VEGFA-Dis-4',
 'VEGFA-Dis-5',
 'VEGFA-Dis-6',
 'VEGFA-Dis-7',
 'VEGFA-Dis-8',
 'VEGFA-OffNo-GDNo-01',
 'VEGFA-OffNo-GDNo-02',
 'VEGFA-OffNo-GDNo-03',
 'VEGFA-OffNo-GDNo-04',
 'VEGFA-OffNo-GDNo-05',
 'VEGFA-OffNo-GDNo-06',
 'VEGFA-OffNo-GDNo-07',
 'VEGFA-OffNo-GDNo-08',
 'VEGFA-OffNo-GDNo-10',
 'VEGFA-OffNo-GDNo-11',
 'VEGFA-OffNo-GDNo-12',
 'VEGFA-OffNo-GDNo-13',
 'VEGFA-OffNo-GDNo-14',
 'VEGFA-OffNo-GDNo-15',
 'VEGFA-OffNo-GDNo-16',
 'VEGFA-OffNot-GDYes-108',
 'VEGFA-OffNot-GDYes-116',
 'VEGFA-OffNot-GDYes-138',
 'VEGFA-OffNot-GDYes-93',
 'VEGFA-OffYes-GDNo-1',
 'VEGFA-OffYes-GDNo-2',
 'VEGFA-OffYes-GDNo-3',
 'VEGFA-OffYes-GDNo-4',
 'VEGFA-OffYes-GDYes-105',
 'VEGFA-OffYes-GDYes-12',
 'VEGFA-OffYes-GDYes-13',
 'VEGFA-OffYes-GDYes-16',
 'VEGFA-OffYes-GDYes-30',
 'VEGFA-OffYes-GDYes-78',
 'VEGFA-off-target-01',
 'VEGFA-off-target-02',
 'VEGFA-off-target-06',
 'VEGFA-off-target-07',
 'VEGFA-off-target-18',
 'VEGFA-off-target-9']
rule all:
    input:
        expand('../TargetSeq_BED_sample_lib/{sample}.bed',sample=SAMPLES),
        expand('../TargetSeq_BED_sample_lib/{sample}.ref.fa',sample=SAMPLES),
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
        for tab in range(1,9):
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
                    bed.write('{chr}\t{tss}\t{tes}\t{flag}\t0\t{strand}\n'.format(chr=chr_,tss=tss,tes=tes,flag=flag,strand=strand))
rule getfasta:
    input:
        '../TargetSeq_BED_sample_lib/{sample}.bed'
    output:
        '../TargetSeq_BED_sample_lib/{sample}.ref.fa'
    shell:
        """
        {BEDTOOLS} getfasta -s -bed {input} -fi {GENOME_HG38} > {output}
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