import os

# ------------------------------------------------------------------->>>>>>>>>>
# Author: Herman ZHAO
# Email: hermanzhaozzzz@gmail.com
# Update
#     date: 2021-07-26 log: format annotation and software path
# ------------------------------------------------------------------->>>>>>>>>>

# ------------------------------------------------------------------->>>>>>>>>>
# SAMPLE INFO
# ------------------------------------------------------------------->>>>>>>>>>
SAMPLES = [
    'ND1-on-target_ND1-L1397N-rep1-0301',
    'ND1-on-target_ND1-L1397N-rep2-0301',
    'ND1-on-target_untreated-rep1-0301',
    'ND1-on-target_untreated-rep2-0301',
    'ND4-on-target_ND4-2000-12-rep1-0301',
    'ND4-on-target_ND4-2000-12-rep2-0301',
    'ND4-on-target_ND4-L1333N-rep1-0301',
    'ND4-on-target_ND4-L1333N-rep2-0301',
    'ND4-on-target_ND4-L1397C-rep1-0301',
    'ND4-on-target_ND4-L1397C-rep2-0301',
    'ND4-on-target_untreated-rep1-0301',
    'ND4-on-target_untreated-rep2-0301',
    'ND5.1-on-target_DddA-NES-ND5.1-rep1-0301',
    'ND5.1-on-target_DddA-NES-ND5.1-rep2-0301',
    'ND5.1-on-target_ND5.1-1333-DddIA-rep1-0301',
    'ND5.1-on-target_ND5.1-1333-DddIA-rep2-0301',
    'ND5.1-on-target_ND5.1-1333-mapKNES-rep1-0301',
    'ND5.1-on-target_ND5.1-1333-rep1-0301',
    'ND5.1-on-target_ND5.1-1333-rep2-0301',
    'ND5.1-on-target_ND5.1-2000-12-rep1-0301',
    'ND5.1-on-target_ND5.1-2000-12-rep2-0301',
    'ND5.1-on-target_TALE-NES-ND5.1-rep1-0301',
    'ND5.1-on-target_TALE-NES-ND5.1-rep2-0301',
    'ND5.1-on-target_UGI-NES-ND5.1-rep1-0301',
    'ND5.1-on-target_UGI-NES-ND5.1-rep2-0301',
    'ND5.1-on-target_WT-ND5.1-rep1-0301',
    'ND5.1-on-target_WT-ND5.1-rep2-0301',
    'ND5.1-on-target_untreated-rep1-0301',
    'ND5.1-on-target_untreated-rep2-0301',
    'ND5.3-on-target_ND5.3-L1397C-rep1-0301',
    'ND5.3-on-target_ND5.3-L1397C-rep2-0301',
    'ND5.3-on-target_untreated-rep1-0301',
    'ND5.3-on-target_untreated-rep2-0301',
    'ND6-on-target_DddA-NES-ND6-rep1',
    'ND6-on-target_DddA-NES-ND6-rep2',
    'ND6-on-target_DddIA-1-0.25-ND6-rep1',
    'ND6-on-target_DddIA-1-0.25-ND6-rep2',
    'ND6-on-target_DddIA-1-0.5-ND6-rep1',
    'ND6-on-target_DddIA-1-0.5-ND6-rep2',
    'ND6-on-target_DddIA-1-1-ND6-rep1',
    'ND6-on-target_DddIA-1-1-ND6-rep2',
    'ND6-on-target_DddIA-1-1.2-ND6-rep1',
    'ND6-on-target_DddIA-1-1.2-ND6-rep2',
    'ND6-on-target_DddIA-1-1.5-ND6-rep1',
    'ND6-on-target_DddIA-1-1.5-ND6-rep2',
    'ND6-on-target_G1309A-ND6-rep1',
    'ND6-on-target_G1309A-ND6-rep2',
    'ND6-on-target_N1308A-ND6-rep1',
    'ND6-on-target_N1308A-ND6-rep2',
    'ND6-on-target_N1367A-ND6-rep1',
    'ND6-on-target_N1367A-ND6-rep2',
    'ND6-on-target_N1368A-ND6-rep1',
    'ND6-on-target_N1368A-ND6-rep2-0301',
    'ND6-on-target_N6-2000-12-rep1-0301',
    'ND6-on-target_N6-2000-12-rep2-0301',
    'ND6-on-target_N6-LTX-12-rep1-0301',
    'ND6-on-target_N6-LTX-12-rep2-0301',
    'ND6-on-target_ND6-1397-DddIA-rep1-0301',
    'ND6-on-target_ND6-1397-DddIA-rep2-0301',
    'ND6-on-target_ND6-1397-hivNES-Q1310A-rep1-0301',
    'ND6-on-target_ND6-1397-hivNES-Q1310A-rep2-0301',
    'ND6-on-target_ND6-1397-hivNES-rep1-0301',
    'ND6-on-target_ND6-1397-hivNES-rep2-0301',
    'ND6-on-target_ND6-1397-rep1-0301',
    'ND6-on-target_ND6-1397-rep2-0301',
    'ND6-on-target_Q1310A-ND6-rep1',
    'ND6-on-target_Q1310A-ND6-rep2',
    'ND6-on-target_TALE-NES-ND6-rep1',
    'ND6-on-target_TALE-NES-ND6-rep2',
    'ND6-on-target_UGI-NES-ND6-rep1',
    'ND6-on-target_UGI-NES-ND6-rep2',
    'ND6-on-target_WT-ND6-rep1',
    'ND6-on-target_WT-ND6-rep2',
    'ND6-on-target_untreated-rep1',
    'ND6-on-target_untreated-rep2',
    'test'
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


# ------------------------------------------------------------------->>>>>>>>>>
# rule all
# ------------------------------------------------------------------->>>>>>>>>>
rule all:
    input:
        expand("../../fastq/{sample}_R1.fastq", sample=SAMPLES),
        expand("../../fastq/{sample}_R2.fastq", sample=SAMPLES),
        expand("../../qc/fastqc/{sample}", sample=SAMPLES),
        "../../qc/multiqc/multiqc_report.html"

# ------------------------------------------------------------------->>>>>>>>>>
# rule fastqc
# ------------------------------------------------------------------->>>>>>>>>>
rule fastqc:
    input:
        "../../fastq/{sample}_R1.fastq",
        "../../fastq/{sample}_R2.fastq"
    output:
        directory("../../qc/fastqc/{sample}")
    log:
        "../../qc/fastqc/{sample}.log"
    shell:
        """
        mkdir -p {output}
        {FASTQC} -o {output} -t {THREADS} \
        {input[0]} {input[1]} > {log} 2>&1
        """

# ------------------------------------------------------------------->>>>>>>>>>
# rule multiqc
# ------------------------------------------------------------------->>>>>>>>>>
rule multiqc:
    input:
        expand("../../qc/fastqc/{sample}",sample=SAMPLES)
    output:
        "../../qc/multiqc/multiqc_report.html"
    log:
        "../../qc/multiqc/multiqc_report.log"
    shell:
        """
        {MULTIQC} ../../qc/fastqc -o ../../qc/multiqc --no-data-dir > {log} 2>&1
        """
