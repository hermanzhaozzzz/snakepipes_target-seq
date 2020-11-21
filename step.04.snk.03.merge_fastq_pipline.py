# _*_ coding: UTF-8 _*_

########################################################################
# Herman ZHAO
# 2020-07-27
# run demultiplex target seq [ APO Mut site ]
########################################################################
# run on abyss
PYTHON = "/home/zhaohuanan/miniconda3/envs/py27/bin/python"




PRIMER_INFO = "./primer_table/sub1.txt"



CUTOFF = ["3","5","10"]

LIBS = ['1341P',
 'DddA-nes',
 'DddIA-1-1-0.5',
 'DddIA-1-1-1',
 'DddIA-1-1-1.5',
 'ND6-WT',
 'TALE-nes',
 'UGI-nes']


rule all:
    input:
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/{lib}_cutoff_{cutoff}.report",lib=LIBS,cutoff=CUTOFF)

rule MergeFastq:
    input:
        "../TargetSeq-{lib}"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/{lib}_cutoff_{cutoff}.report"
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/{lib}_cutoff_{cutoff}.log"
    params:
        out_dir = "../TargetSeq-{lib}/cutoff_{cutoff}",
        cutoff = "{cutoff}"
    shell:
#         """
#         srun -T 24 \
        """
        {PYTHON} ./program/target_seq_merge_fq_V03.py \
        -i {input[0]} \
        -p {PRIMER_INFO} \
        -o {params.out_dir} \
        -r {output} \
        --MinMergeReadNumCutoff {params.cutoff}  > {log} 2>&1
        """

