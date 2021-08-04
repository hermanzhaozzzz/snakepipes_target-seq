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
LIB = ['ABE-HEK4-rep1',
 'ABE-HEK4-rep2',
 'ABE8e-HEK4-rep1',
 'ABE8e-HEK4-rep2',
 'ACBE-HEK4-rep1',
 'ACBE-HEK4-rep2',
 'CBE-HEK4-rep1',
 'CBE-HEK4-rep2']
REGION = ['HEK4-ABEABE8EACBECBEOUT-Detect1',
 'HEK4-ABEABE8EACBECBEOUT-Digenome34',
 'HEK4-ABEABE8EACBECBEOUT-N2',
 'HEK4-ABEABE8EACBEcbeOUT-C1',
 'HEK4-ABEABE8EACBEcbeout-B2',
 'HEK4-ABEABE8Eacbecbeout-D1',
 'HEK4-ABEABE8Eacbecbeout-D3',
 'HEK4-ABEABE8Eacbecbeout-D4',
 'HEK4-abeABE8EACBECBEout-A1',
 'HEK4-abeABE8EACBECBEout-A2',
 'HEK4-abeABE8EACBECBEout-A6',
 'HEK4-abeabe8eACBECBEout-A5',
 'HEK4-abeabe8eacbeCBEout-A8']

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