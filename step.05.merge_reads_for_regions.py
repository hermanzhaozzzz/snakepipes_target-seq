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
LIBS = ['untreated-rep1', 'untreated-rep2']


CUTOFF = [
    "0","3","5","10"
         ]


# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
THREADS = '20'

# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
PRIMER_INFO = "./primer_table/primer.txt"

# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# polaris
PYTHON = "/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/miniconda3/envs/snakepipes_py27/bin/python"
# abyss
# PYTHON = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq/bin/python"


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
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
        if (({wildcards.cutoff}==0)); then 
        echo "cutoff == 0!!" 
        mkdir -p {input}/cutoff_0/merge.fastq
        rsync -Ptrv {input}/demultiplex.fastq/* {input}/cutoff_0/merge.fastq
        cd {input}/demultiplex.fastq
        zsh ../../snakepipes_target-seq/cutoff_0_rename.sh
        cd -
        touch {output}
        else
        {PYTHON} ./program/target_seq_merge_fq_V03.py \
        -i {input[0]} \
        -p {PRIMER_INFO} \
        -o {params.out_dir} \
        -r {output} \
        --MinMergeReadNumCutoff {params.cutoff}  > {log} 2>&1
        fi
        """
