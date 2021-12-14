import os

# ------------------------------------------------------------------->>>>>>>>>>
# Author: Herman ZHAO
# Email: hermanzhaozzzz@gmail.com
# Update
#     date: 2021-08-03 log: format annotation and software path
# ------------------------------------------------------------------->>>>>>>>>>

# ------------------------------------------------------------------->>>>>>>>>>
# SAMPLE INFO
# ------------------------------------------------------------------->>>>>>>>>>
LIB = ['untreated-rep1', 'untreated-rep2']

REGION = ['ND4-new-only1',
 'ND4-new-only10',
 'ND4-new-only11',
 'ND4-new-only12',
 'ND4-new-only13',
 'ND4-new-only2',
 'ND4-new-only3',
 'ND4-new-only4',
 'ND4-new-only5',
 'ND4-new-only6',
 'ND4-new-only7',
 'ND4-new-only8',
 'ND4-new-only9',
 'ND5.1-new-only1',
 'ND5.1-new-only10',
 'ND5.1-new-only11',
 'ND5.1-new-only12',
 'ND5.1-new-only13',
 'ND5.1-new-only14',
 'ND5.1-new-only15',
 'ND5.1-new-only16',
 'ND5.1-new-only17',
 'ND5.1-new-only18',
 'ND5.1-new-only19',
 'ND5.1-new-only2',
 'ND5.1-new-only20',
 'ND5.1-new-only21',
 'ND5.1-new-only22',
 'ND5.1-new-only23',
 'ND5.1-new-only24',
 'ND5.1-new-only25',
 'ND5.1-new-only3',
 'ND5.1-new-only4',
 'ND5.1-new-only5',
 'ND5.1-new-only6',
 'ND5.1-new-only7',
 'ND5.1-new-only8',
 'ND5.1-new-only9',
 'ND6-new-only1',
 'ND6-new-only10',
 'ND6-new-only11',
 'ND6-new-only12',
 'ND6-new-only13',
 'ND6-new-only14',
 'ND6-new-only15',
 'ND6-new-only16',
 'ND6-new-only17',
 'ND6-new-only2',
 'ND6-new-only3',
 'ND6-new-only4',
 'ND6-new-only5',
 'ND6-new-only6',
 'ND6-new-only7',
 'ND6-new-only8',
 'ND6-new-only9',
 'share-new-1',
 'share-new-10',
 'share-new-11',
 'share-new-12',
 'share-new-13',
 'share-new-14',
 'share-new-15',
 'share-new-16',
 'share-new-17',
 'share-new-18',
 'share-new-19',
 'share-new-2',
 'share-new-20',
 'share-new-21',
 'share-new-22',
 'share-new-23',
 'share-new-24',
 'share-new-25',
 'share-new-26',
 'share-new-27',
 'share-new-28',
 'share-new-29',
 'share-new-3',
 'share-new-30',
 'share-new-31',
 'share-new-32',
 'share-new-4',
 'share-new-5',
 'share-new-6',
 'share-new-7',
 'share-new-8',
 'share-new-9']

# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
THREADS = '20'

# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# none

# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# get the application path
with os.popen("which fastqc") as path:
    FASTQC = path.read().strip()
    print('PATH fastqc:', FASTQC)
with os.popen("which multiqc") as path:
    MULTIQC = path.read().strip()
    print('PATH multiqc:', MULTIQC)

if "" in [FASTQC, MULTIQC]:
    raise ValueError("The necessary software path is missing!!!")
else:
    print("pass software check...\n")







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
    log:
        "../qc_quantify-all/fastqc/TargetSeq-{lib}/{region}.log"
    shell:
        """
        mkdir -p {output}
        {FASTQC} -o {output} -t {THREADS} \
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