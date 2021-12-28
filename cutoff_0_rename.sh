for i in `ls ../demultiplex.fastq | grep fastq`                                               
mv ../cutoff_0/merge.fastq/$i ../cutoff_0/merge.fastq/${${i}/_demultiplex_/_merge_barcode_};


echo "rename done!"