# ------------------------------------------------------------------->>>>>>>>>>
# 需要编辑的变量
# ------------------------------------------------------------------->>>>>>>>>>
# polaris
CONDA_ENV=/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/miniconda3/envs/snakepipes_py27
# abyss
# CONDA_ENV=/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq

PRIMER_INFO=./primer_table/primer.txt
# ------------------------------------------------------------------->>>>>>>>>>
# 脚本内容
# ------------------------------------------------------------------->>>>>>>>>>
PYTHON2=${CONDA_ENV}/bin/python

for case_name in untreated-rep1 untreated-rep2
do    
    fq1=../fastq/${case_name}_R1.fastq.gz
    fq2=../fastq/${case_name}_R2.fastq.gz
    out_dir=../IntactTargetSeq-${case_name}
    pkurun-cnnl 1 20 $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py -1 $fq1 -2 $fq2 -p $PRIMER_INFO -o $out_dir
    sleep 1
done











