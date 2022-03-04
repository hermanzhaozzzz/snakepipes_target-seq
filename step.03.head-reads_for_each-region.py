import glob
# 清除文件
# rm ../Int*/*/*.fastq.gz
# rm ../Target*/*/*.fastq.gz




SAMPLE = [
    'ND1-Det-Q5U-rep1',
    'ND1-Det-Q5U-rep2',
    'ND4-Det-Q5U-rep1',
    'ND4-Det-Q5U-rep2',
    'ND4-L1333N-Det-Q5U-rep1',
    'ND4-L1333N-Det-Q5U-rep2',
    'ND4-L1397C-Det-Q5U-rep1',
    'ND4-L1397C-Det-Q5U-rep2',
    'ND5-1-Det-Q5U-rep1',
    'ND5-1-Det-Q5U-rep2',
    'ND5-3-L1397C-Det-Q5U-rep1',
    'ND5-3-L1397C-Det-Q5U-rep2',
    'ND6-Det-Q5U-rep1',
    'ND6-Det-Q5U-rep2',
    'ND6-LTX-12-Q5U-rep1',
    'ND6-LTX-12-Q5U-rep2'
]

FQS = []

for lib in SAMPLE:
    ls_this_folder = glob.glob(f'../IntactTargetSeq-{lib}/demultiplex.fastq/*.fastq')
    ls_fq = [f'{lib}/demultiplex.fastq/{fastq.split("demultiplex.fastq/")[1]}' for fastq in ls_this_folder]
    FQS.extend(ls_fq)

print(FQS)

# 取多少行？
# 6 000 000 lines
# 1.5M reads
# 3M PE reads
LINES = 10000000
THREADS = 20

rule all:
    input:
        expand("../TargetSeq-{fq}.gz", fq=FQS)
rule head:
    input:
        "../IntactTargetSeq-{fq}"
    output:
        temp("../TargetSeq-{fq}")
    shell:
        """
        head -{LINES} {input} > {output}
        pigz --keep -p {THREADS} {input}
        """
rule gzip_output:
    input:
        "../TargetSeq-{fq}"
    output:
        "../TargetSeq-{fq}.gz"
    shell:
        """
        pigz -p {THREADS} {input}
        """

