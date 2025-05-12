PRIMER_INFO=./primer_table/primer.txt
# ------------------------------------------------------------------->>>>>>>>>>
# 脚本内容
# ------------------------------------------------------------------->>>>>>>>>>

for case_name in Q5U-56-rep1 Q5U-56-rep2 Q5U-mock-rep1 Q5U-mock-rep2
do    
    fq1=../fastq/${case_name}_R1.fastq.gz
    fq2=../fastq/${case_name}_R2.fastq.gz
    out_dir=../IntactTargetSeq-${case_name}

    
    python3 ./program/demultiplex.py \
        --fastq_read1 $fq1 \
        --fastq_read2 $fq2 \
        --primer_info_file $PRIMER_INFO \
        --output_dir $out_dir
    sleep 1
done











