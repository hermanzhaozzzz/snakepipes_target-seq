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
LIB = ['Q5U-56-rep1', 'Q5U-56-rep2', 'Q5U-mock-rep1', 'Q5U-mock-rep2']



REGION = [
    'ND5.1-new-only11',
    'ND5.1-new-only14',
    'ND5.1-new-only17',
    'ND6-new-only10',
    'ND6-new-only14',
    'ND6-new-only16',
    'ND6-new-only17',
    'ND6-new-only2',
    'ND6-new-only4',
    'share-new-12',
    'share-new-30',
    'share-new-6'
]

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
        expand("../IntactTargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R1.fastq.gz",lib=LIB, region=REGION),
        expand("../IntactTargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R2.fastq.gz",lib=LIB, region=REGION),
        expand("../qc_quantify-all/fastqc/TargetSeq-{lib}_{region}",lib=LIB, region=REGION),
        "../qc_quantify-all/multiqc/multiqc_report.html"
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule fastqc
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule cp:
    input: 
        "../IntactTargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R1.fastq.gz",
        "../IntactTargetSeq-{lib}/demultiplex.fastq/{region}_demultiplex_R2.fastq.gz"
    output: 
        temp("../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R1.fastq.gz"),
        temp("../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R2.fastq.gz"),
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """
rule fastqc: 
    input: 
        "../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R1.fastq.gz",
        "../qc_quantify-all/fastq/TargetSeq-{lib}_{region}_demultiplex_R2.fastq.gz",
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