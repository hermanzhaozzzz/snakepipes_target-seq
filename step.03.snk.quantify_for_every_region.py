# _*_ coding: UTF-8 _*_

########################################################################
# ZHAO Huanan
# 2020-05-04
# fastqc + multiqc pipeline
######################################################################## 
# run on abyss

# --------------------------------------------------------------->>>>>>>
# pipeline
# --------------------------------------------------------------->>>>>>>
# 1. fastqc
# 2. multiqc
# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
FASTQC = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/fastqc"
MULTIQC = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/multiqc"





# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
LIB = ['N1-veri-2', 'N4-1333-2', 'N4-1397R-2', 'N4-veri-2', 'N5-3-veri-2']
REGION = ['ND516-share-10',
 'ND516-share-11',
 'ND516-share-12',
 'ND516-share-13',
 'ND516-share-2',
 'ND516-share-3',
 'ND516-share-4',
 'ND516-share-5',
 'ND516-share-6',
 'ND516-share-7',
 'ND516-share-8',
 'ND516-share-9']



# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    # 有点特殊，输出和输入文件名只能有R1.fastq.gz和_fastqc.html的区别，前面不能有区别！要改就得全改！
    input:
        expand("../TargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R1.fastq",lib=LIB, region=REGION),
        expand("../TargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R2.fastq",lib=LIB, region=REGION),
        expand("../qc_quantify-all/fastqc/TargetSeq-{lib}_{region}",lib=LIB, region=REGION),
        "../qc_quantify-all/multiqc/multiqc_report.html"
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule fastqc
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule cp:
    input: 
        "../TargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R1.fastq",
        "../TargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R2.fastq"
    output: 
        temp("../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R1.fastq"),
        temp("../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R2.fastq"),
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """
rule fastqc: 
    input: 
        "../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R1.fastq",
        "../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R2.fastq",
    output: 
        directory("../qc_quantify-all/fastqc/TargetSeq-{lib}_{region}")
                # the suffix _fastqc.zip is necessary the same with SAMPLES raw name suffix, for multiqc to find the file. 
                # If not using multiqc, you are free to choose an arbitrary filename
    params: 
        threads = "24",
    log:
        "../qc_quantify-all/fastqc/TargetSeq-{lib}/{region}.log"
    shell:
        """
        mkdir -p {output}
        {FASTQC} -o {output} -t {params.threads} \
        {input[0]} {input[1]} 2>{log}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule multiqc
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule multiqc: 
    input: 
        expand("../qc_quantify-all/fastqc/TargetSeq-{lib}_{region}",lib=LIB, region=REGION)
                # the suffix XXX_fastqc.zip is necessary the same with SAMPLES raw name suffix, for multiqc to find the file. 
    output:
        "../qc_quantify-all/multiqc/multiqc_report.html"
    log:
        "../qc_quantify-all/multiqc.log"
    shell: 
        """
        {MULTIQC} ../qc_quantify-all/fastqc -o ../qc_quantify-all/multiqc --no-data-dir 2>{log}
        """