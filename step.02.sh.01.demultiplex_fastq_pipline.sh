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
PYTHON2=/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq/bin/python

# ########################################################################

for case_name in N4-veri-2 N1-veri-2 N5-3-veri-2 N4-1333-2 N4-1397R-2
do
    PRIMER_INFO=./primer_table/checking.txt
    
    
    fq1=../fastq/${case_name}_combined_R1.fastq.gz
    fq2=../fastq/${case_name}_combined_R2.fastq.gz
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