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
PRIMER_INFO=./primer_table/8.txt

for case_name in M8-B-1 M8-B-2 M9-Y-1 M9-Y-2
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
