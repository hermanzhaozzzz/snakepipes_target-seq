# _*_ coding: UTF-8 _*_

########################################################################
# Herman ZHAO
# 2020-07-27
# run demultiplex target seq [ APO Mut site ]
########################################################################
# run on abyss
PYTHON = "/home/zhaohuanan/miniconda3/envs/py27/bin/python"




PRIMER_INFO = "./primer_table/1.txt"



CUTOFF = ["3"]

LIBS = [
    'B-1',
    'B-2',
    'M1-1',
    'M1-2',
    'M2-1',
    'M2-2', 
    'M3-1',
    'M3-2',
    'M4-1',
    'M4-2',
    'M5-1',
    'M5-2',
    'M6-1',
    'M6-2',
    'M7-1',
    'M7-2',
    'S334-1',
    'S334-2',
    'Y-1',
    'Y-2'
]


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
        """
        srun -T 24 \
        {PYTHON} ./program/target_seq_merge_fq_V03.py \
        -i {input[0]} \
        -p {PRIMER_INFO} \
        -o {params.out_dir} \
        -r {output} \
        --MinMergeReadNumCutoff {params.cutoff}  > {log} 2>&1
        """
