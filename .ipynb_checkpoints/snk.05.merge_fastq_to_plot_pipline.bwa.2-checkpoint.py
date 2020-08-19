# _*_ coding: UTF-8 _*_

########################################################################
# MENG Howard
# 2020-07-07
# run demultiplex target seq [ APO Mut site ]
########################################################################
# run on abyss
# /home/menghaowei/menghw_HD/BE_project/20.target_seq_all/16.APOmut_screen.20200707
CUTADAPT = "/home/zhaohuanan/miniconda3/envs/cutadapt/bin/cutadapt"
BWA = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bwa"
SAMTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_cutadapt-STARmapping-FPKM-sortBAM/bin/samtools"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/bedtools"
SAMCLIP = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/samclip"
PYTHON = "/home/zhaohuanan/miniconda3/envs/snakepipes_target-seq-from-table-to-plot/bin/python"





CUTOFF = ["3"]

LIBS = [
    'EH-1',
    'EH-2'
]

SAMPLES = ['EMX1-notOFF-01',
 'EMX1-notOFF-02',
 'EMX1-notOFF-03',
 'HEK4-Guideseq-3',
 'HEK4-Guideseq-4',
 'HEK4-Guideseq-6']

READ_IDX = ["1","2"]

defult_sgRNA_dict_for_plot = {
    'VEGFA': "GACCCCCTCCACCCCGCCTCCGG", 
    'EMX1': "GAGTCCGAGCAGAAGAAGAAGGG", 
    'HEK3': "GGCCCAGACTGAGCACGTGATGG", 
    "HEK4":"GGCACTGCGGCTGGAGGTGGGGG", 
    "RNF2":"GTCATCTTAGTCATTACCTGAGG"
}

rule all:
    input:
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R{read_idx}_cutadapt.fq.gz",lib=LIBS,sample=SAMPLES,read_idx=READ_IDX,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF),
        expand("../all_plot/cutoff_{cutoff}.ext50/TargetSeq-{lib}_{sample}_cutoff_{cutoff}_indel.ext50.pdf",lib=LIBS,sample=SAMPLES,cutoff=CUTOFF)
        

# wildcard_constraints:
#     pass


rule check_file:
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq",
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R2.fastq"
    run:
        try:
            open(output[0],'r')
        except IOError:
            open(output[0],'w')
        try:
            open(output[1],'r')
        except IOError:
            open(output[1],'w')
            
rule cutadapt:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R1.fastq",
        "../TargetSeq-{lib}/cutoff_{cutoff}/merge.fastq/{sample}_merge_barcode_R2.fastq"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz"
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_cutadapt.log"
    shell:
        #"""
        #srun -T 24 \
        """
        {CUTADAPT} -j 24 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
        """  

rule bwa_mapping:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R1_cutadapt.fq.gz",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_R2_cutadapt.fq.gz",
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam"
    params:
        ref = '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
    log:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.log"
    shell:
        """
        {BWA} mem {params.ref} {input[0]} {input[1]} -t 24 -M > {output} 2>{log}
        """
rule sam_to_bam:
    input:
        sam = "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.sam",
        ref = '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    shell:
        """
        {SAMTOOLS} view -h -f 1 -F 268 {input.sam} \
        | {SAMCLIP} --ref {input.ref} --max 3 --progress 0 \
        | awk 'function abs(v) {{return v < 0 ? -v : v}} $1~"@" || ($7 == "=" && abs($9) <= 600 ) {{print $0}}' \
        | {SAMTOOLS} view -hb > {output}
        """


rule samtools_sort_by_position:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa.bam"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    shell:
        "{SAMTOOLS} sort -O BAM -o {output} -T {output}.temp -@ 6 -m 2G {input}"


rule samtools_index:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai"
    shell:
        "{SAMTOOLS} index -@ 6 {input} {output}"


rule spike_in_mpileup:
    input:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam",
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bam.bai"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup"
    params:
        ref = '../TargetSeq_BED_sample_lib/{sample}.ref.upper.fa'
    shell:
        "{SAMTOOLS} mpileup {input[0]} --reference {params.ref} --max-depth 10000000 -q 20 -Q 20 > {output} " 


rule parse_mpileup:
    input: 
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.mpileup"
    output:
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat"
    shell:
        "{PYTHON} ./program/parse-mpileup-V04.py -i {input} -o {output} -n 0"


rule bmat_plot:
    input: 
        "../TargetSeq-{lib}/cutoff_{cutoff}/mapping/{sample}_bwa_sort.bmat"
    output:
        "../all_plot/cutoff_{cutoff}.ext50/TargetSeq-{lib}_{sample}_cutoff_{cutoff}_indel.ext50.pdf"
    params:
        sgRNA_seq = lambda wildcards, input: defult_sgRNA_dict_for_plot[input[0].split("/")[4].split("-")[0]]
    shell:
        "{PYTHON} ./program/plot-targetseq-bmat-V02.py -i {input} -o {output} --region_extend_length 50 --sgRNA {params.sgRNA_seq}"
######################################################################## 
# separate plot
######################################################################## 
# import os

# defult_sgRNA_dict = {
#     'VEGFA': "GACCCCCTCCACCCCGCCTCCGG", 
#     'EMX1': "GAGTCCGAGCAGAAGAAGAAGGG", 
#     'HEK3': "GGCCCAGACTGAGCACGTGATGG", 
#     "HEK4":"GGCACTGCGGCTGGAGGTGGGGG", 
#     "RNF2":"GTCATCTTAGTCATTACCTGAGG"
# }   

# for lib in LIBS:
#     cmd_fmt = 'plot-targetseq-bmat-V02.py  -i {input} -o {output} --region_extend_length 80 --sgRNA {sgRNA_seq} &'

#     for sample in SAMPLES:
#         bmat_file = "293T-TargetSeq-{lib_name}/cutoff_3/mapping/{sample_name}_bwa_sort.bmat".format(lib_name=lib, sample_name=sample)
#         sgRNA_seq = defult_sgRNA_dict.get(sample.split("-")[0])

#         # make output
#         out_filename = "all_plot/cutoff_3/293T-TargetSeq-{lib_name}_{sample_name}_cutoff_3_indel.ext80.E_sgRNA.pdf".format(lib_name=lib, sample_name=sample)
    
#         # make cmd
#         cmd = cmd_fmt.format(
#             input = bmat_file,
#             output = out_filename,
#             sgRNA_seq = sgRNA_seq)

#         print(cmd)
#         os.system(cmd)


######################################################################## 
# merge plot 
######################################################################## 
# import os

# defult_sgRNA_dict = {
#     'VEGFA': "GACCCCCTCCACCCCGCCTCCGG", 
#     'EMX1': "GAGTCCGAGCAGAAGAAGAAGGG", 
#     'HEK3': "GGCCCAGACTGAGCACGTGATGG", 
#     "HEK4":"GGCACTGCGGCTGGAGGTGGGGG", 
#     "RNF2":"GTCATCTTAGTCATTACCTGAGG"
# }   

# # mkdir 
# os.system("mkdir -p all_plot/all_comparison/")

# # tools
# plot_code = "~/menghw_HD/BE_project/20.target_seq_all/tools/target_seq_plot_rep_v5.py"

# # run
# for sample in SAMPLES:
#     cmd_fmt = 'python2 {plot_code} -i {input_bmat_file} -n {lib_name} -o {out_plot} -s {sgRNA_seq} --indel False --region 25 --format pdf &'
#     bmat_file_list = [] 

#     for lib in LIBS:
#         bmat_file = "293T-TargetSeq-{lib_name}/cutoff_3/mapping/{sample_name}_bwa_sort.bmat".format(lib_name=lib, sample_name=sample)
#         bmat_file_list.append(bmat_file)

#     # make cmd 
#     input_bmat_file = ",".join(bmat_file_list)
#     lib_name = ",".join(LIBS)
#     out_filename = "all_plot/all_comparison/{sample_name}_mutant_cutoff3.pdf".format(sample_name=sample)

#     sgRNA_seq = defult_sgRNA_dict.get(sample.split("-")[0])

#     if sgRNA_seq == None:
#         sgRNA_seq = "GGCACTGCGGCTGGAGGTGGGGG"

#     cmd = cmd_fmt.format(
#         plot_code = plot_code,
#         input_bmat_file = input_bmat_file,
#         lib_name = lib_name,
#         out_plot = out_filename,
#         sgRNA_seq = sgRNA_seq)

#     print(cmd)
    # os.system(cmd)



