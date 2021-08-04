# ------------------------------------------------------------------->>>>>>>>>>
# 需要编辑的变量
# ------------------------------------------------------------------->>>>>>>>>>
# polaris
CONDA_ENV=/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/miniconda3/envs/snakepipes_py27
# abyss
# CONDA_ENV=/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq

PRIMER_INFO=./primer_table/ABE.txt
# ------------------------------------------------------------------->>>>>>>>>>
# 脚本内容
# ------------------------------------------------------------------->>>>>>>>>>
PYTHON2=${CONDA_ENV}/bin/python

for case_name in ABE-HEK4-rep1 ABE-HEK4-rep2 ABE8e-HEK4-rep1 ABE8e-HEK4-rep2 ACBE-HEK4-rep1 ACBE-HEK4-rep2 CBE-HEK4-rep1 CBE-HEK4-rep2
do    
    fq1=../fastq/${case_name}_combined_R1.fastq.gz
    fq2=../fastq/${case_name}_combined_R2.fastq.gz
    out_dir=../TargetSeq-${case_name}
    pkurun-cns 1 20 $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py -1 $fq1 -2 $fq2 -p $PRIMER_INFO -o $out_dir
    sleep 1
done











