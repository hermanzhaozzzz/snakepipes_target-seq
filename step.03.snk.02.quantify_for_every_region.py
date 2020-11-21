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

# conda install -c bioconda fastqc=0.11.8
FASTQC = "~/miniconda3/envs/snakepipes_fastqc-multiqc/bin/fastqc"

# conda install -c bioconda multiqc=1.7
MULTIQC = "~/miniconda3/envs/snakepipes_fastqc-multiqc/bin/multiqc"





# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
LIB = ['1341P',
 'DddA-nes',
 'DddIA-1-1-0.5',
 'DddIA-1-1-1',
 'DddIA-1-1-1.5',
 'ND6-WT',
 'TALE-nes',
 'UGI-nes']
REGION = ['ND516-share-4',
 'ND516-share-5',
 'ND516-share-9',
 'ND6-only-11',
 'ND6-only-2',
 'ND6-only-6']




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
        "../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R1.fastq",
        "../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R2.fastq",
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