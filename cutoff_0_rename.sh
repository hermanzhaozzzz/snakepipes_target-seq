for i in `ls ../demultiplex.fastq | grep fastq`;
do
    mv ../cutoff_0/merge.fastq/$i ../cutoff_0/merge.fastq/${${i}/_demultiplex_/_merge_barcode_};
done

echo "rename done!"