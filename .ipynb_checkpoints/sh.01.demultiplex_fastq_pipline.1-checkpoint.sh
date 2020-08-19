# _*_ coding: UTF-8 _*_

########################################################################
# MENG Howard
# 2020-07-07
# run demultiplex target seq [ APO Mut site ]
########################################################################
# run on abyss
# /home/menghaowei/menghw_HD/BE_project/20.target_seq_all/16.APOmut_screen.20200707



# 先salloc节点，再跑
PYTHON2=/home/zhaohuanan/miniconda3/envs/py27/bin/python
PRIMER_INFO=./primer_table/1.txt

for case_name in B-1 B-2 M1-1 M1-2 M2-1 M2-2 M3-1 M3-2 M4-1 M4-2 M5-1 M5-2 M6-1 M6-2 M7-1 M7-2 S334-1 S334-2 Y-1 Y-2
do

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
echo "all done"
