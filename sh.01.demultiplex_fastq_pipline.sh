# _*_ coding: UTF-8 _*_

########################################################################
# ZHAO Huanan
# 2020-10-10
# run demultiplex target seq [ APO Mut site ]
########################################################################
# run on abyss
# /gpfs/user/zhaohuanan/3.project/13.2020-10_1010_TargetSeq-GBE-DdCBE/snakepipes_target-seq-from-table-to-plot
########################################################################



# 先salloc节点，再跑
PYTHON2=/home/zhaohuanan/miniconda3/envs/py27/bin/python

########################################################################

for case_name in 5 6 7
do
    PRIMER_INFO=./primer_table/5-6-7.txt
    
    
    fq1=../fix.fastq/TargetSeq-${case_name}_R1.fastq.gz
    fq2=../fix.fastq/TargetSeq-${case_name}_R2.fastq.gz
    out_dir=../TargetSeq-${case_name}

    mkdir -p $out_dir
    log=../TargetSeq-${case_name}/TargetSeq-${case_name}_site_split_FASTQ.log

    $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py \
    -1 $fq1 -2 $fq2 \
    -p $PRIMER_INFO \
    -o $out_dir > $log 2>&1 &
    echo $case_name "done"
done
########################################################################
for case_name in 8 9 10
do
    PRIMER_INFO=./primer_table/8-9-10.txt
    
    
    fq1=../fix.fastq/TargetSeq-${case_name}_R1.fastq.gz
    fq2=../fix.fastq/TargetSeq-${case_name}_R2.fastq.gz
    out_dir=../TargetSeq-${case_name}

    mkdir -p $out_dir
    log=../TargetSeq-${case_name}/TargetSeq-${case_name}_site_split_FASTQ.log

    $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py \
    -1 $fq1 -2 $fq2 \
    -p $PRIMER_INFO \
    -o $out_dir > $log 2>&1 &
    echo $case_name "done"
done
########################################################################
for case_name in 11 12 13
do
    PRIMER_INFO=./primer_table/11-12-13.txt
    
    
    fq1=../fix.fastq/TargetSeq-${case_name}_R1.fastq.gz
    fq2=../fix.fastq/TargetSeq-${case_name}_R2.fastq.gz
    out_dir=../TargetSeq-${case_name}

    mkdir -p $out_dir
    log=../TargetSeq-${case_name}/TargetSeq-${case_name}_site_split_FASTQ.log

    $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py \
    -1 $fq1 -2 $fq2 \
    -p $PRIMER_INFO \
    -o $out_dir > $log 2>&1 &
    echo $case_name "done"
done
########################################################################
for case_name in 14 15
do
    PRIMER_INFO=./primer_table/14-15.txt
    
    
    fq1=../fix.fastq/TargetSeq-${case_name}_R1.fastq.gz
    fq2=../fix.fastq/TargetSeq-${case_name}_R2.fastq.gz
    out_dir=../TargetSeq-${case_name}

    mkdir -p $out_dir
    log=../TargetSeq-${case_name}/TargetSeq-${case_name}_site_split_FASTQ.log

    $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py \
    -1 $fq1 -2 $fq2 \
    -p $PRIMER_INFO \
    -o $out_dir > $log 2>&1 &
    echo $case_name "done"
done
########################################################################
for case_name in duigou jiantou
do
    PRIMER_INFO=./primer_table/duigou-jiantou.txt
    
    
    fq1=../fix.fastq/TargetSeq-${case_name}_R1.fastq.gz
    fq2=../fix.fastq/TargetSeq-${case_name}_R2.fastq.gz
    out_dir=../TargetSeq-${case_name}

    mkdir -p $out_dir
    log=../TargetSeq-${case_name}/TargetSeq-${case_name}_site_split_FASTQ.log

    $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py \
    -1 $fq1 -2 $fq2 \
    -p $PRIMER_INFO \
    -o $out_dir > $log 2>&1 &
    echo $case_name "done"
done
########################################################################
for case_name in quanquan sanjiao
do
    PRIMER_INFO=./primer_table/quanquan-sanjiao.txt
    
    
    fq1=../fix.fastq/TargetSeq-${case_name}_R1.fastq.gz
    fq2=../fix.fastq/TargetSeq-${case_name}_R2.fastq.gz
    out_dir=../TargetSeq-${case_name}

    mkdir -p $out_dir
    log=../TargetSeq-${case_name}/TargetSeq-${case_name}_site_split_FASTQ.log

    $PYTHON2 ./program/target_seq_demultiplex_fastq_V01.py \
    -1 $fq1 -2 $fq2 \
    -p $PRIMER_INFO \
    -o $out_dir > $log 2>&1 &
    echo $case_name "done"
done
########################################################################


echo "all done"