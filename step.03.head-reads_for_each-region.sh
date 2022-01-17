#!/bin/bash
for case_name in Q5U-56-rep1 Q5U-56-rep2 Q5U-mock-rep1 Q5U-mock-rep2
do    
    mkdir -p ../TargetSeq-${case_name}/demultiplex.fastq;
    rm -rf ../TargetSeq-${case_name}/demultiplex.fastq/*;
    
    for fastq in `ls ../IntactTargetSeq-${case_name}/* | grep fastq`;
    do
        head -10000000 ../IntactTargetSeq-${case_name}/demultiplex.fastq/$fastq > ../TargetSeq-${case_name}/demultiplex.fastq/$fastq;
        pigz -p 24 ../IntactTargetSeq-${case_name}/demultiplex.fastq/$fastq;
        pigz -p 24 ../TargetSeq-${case_name}/demultiplex.fastq/$fastq;
        echo ${case_name} ${fastq} "head -10000000 Done"
    done
    
done
echo "All Done"

# 6 000 000 lines
# 1.5M reads
# 3M PE reads
