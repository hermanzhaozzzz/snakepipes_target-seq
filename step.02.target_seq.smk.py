# snakemake -p -j 50 -s step.02.demultiplex_fastq.smk --rerun-incomplete --cores 8
configfile: "config.yaml"

# 从配置中获取参数
cases = config["cases"]
regions = config["regions"]
primer_info = config["primer_info"]
input_dir = config["input_dir"]
output_base = config["output_base"]
mismatch = config["mismatch"]
max_reads = config["max_reads"]
merge_cutoff = config["merge_cutoff"]
genome_hg38 = config["genome_hg38"]
cutadapt_r1 = config["cutadapt_r1"]
cutadapt_r2 = config["cutadapt_r2"]
mut_direction = config["mut_direction"]
rule all:
    input:
        # Demultiplex输出文件（按case+region展开）
        expand("{output_base}/demultiplex/{case}/{region}_R1.fastq.gz", output_base=output_base, case=cases, region=regions),
        expand("{output_base}/demultiplex/{case}/{region}_R2.fastq.gz", output_base=output_base, case=cases, region=regions),
        # FastQC报告
        expand("{output_base}/qc/fastqc/{case}/{region}_R1_fastqc.html", output_base=output_base, case=cases, region=regions),
        expand("{output_base}/qc/fastqc/{case}/{region}_R2_fastqc.html", output_base=output_base, case=cases, region=regions),
        # # MultiQC报告
        expand("{output_base}/qc/multiqc_report.html", output_base=output_base),
        # 合并后的fastq文件
        expand("{output_base}/merged/{case}/{region}_cutoff_{cutoff}_R1.fastq.gz", output_base=output_base, case=cases, region=regions, cutoff=merge_cutoff),
        expand("{output_base}/merged/{case}/{region}_cutoff_{cutoff}_R2.fastq.gz", output_base=output_base, case=cases, region=regions, cutoff=merge_cutoff),
        expand("{output_base}/merged/{case}/{region}_cutoff_{cutoff}.report", output_base=output_base, case=cases, region=regions, cutoff=merge_cutoff),
        # # relinked文件
        expand("{output_base}/merged_groupped_by_region/{region}/{region}_{case}_cutoff_{cutoff}_R1.fastq.gz", output_base=output_base, case=cases, region=regions, cutoff=merge_cutoff),
        expand("{output_base}/merged_groupped_by_region/{region}/{region}_{case}_cutoff_{cutoff}_R2.fastq.gz", output_base=output_base, case=cases, region=regions, cutoff=merge_cutoff),
        # # 参考基因组索引文件
        expand('{output_base}/refs/{region}.ref.upper.fa', output_base=output_base,region=regions),
        expand('{output_base}/refs/{region}.ref.upper.fa.fai', output_base=output_base,region=regions),
        expand('{output_base}/refs/{region}.ref.upper.fa.amb', output_base=output_base,region=regions),
        expand('{output_base}/refs/{region}.ref.upper.fa.ann', output_base=output_base,region=regions),
        expand('{output_base}/refs/{region}.ref.upper.fa.bwt', output_base=output_base,region=regions),
        expand('{output_base}/refs/{region}.ref.upper.fa.pac', output_base=output_base,region=regions),
        expand('{output_base}/refs/{region}.ref.upper.fa.sa', output_base=output_base,region=regions),
        expand('{output_base}/refs/{region}.focus.upper.fa.seq', output_base=output_base,region=regions),
        # # plot region 文件
        expand("{output_base}/plots/single/{region}.{case}.cutoff_{cutoff}.base_mut.tsv", output_base=output_base, case=cases, region=regions, cutoff=merge_cutoff),
        expand("{output_base}/plots/single/{region}.{case}.cutoff_{cutoff}.base_mut.plot.pdf", output_base=output_base, case=cases, region=regions, cutoff=merge_cutoff),
        # # plot multi 文件
        expand("{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.heatmap.pdf", output_base=output_base, region=regions, cutoff=merge_cutoff),
        expand("{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.heatmap.csv", output_base=output_base, region=regions, cutoff=merge_cutoff),
        expand("{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.count_ratio.pdf", output_base=output_base, region=regions, cutoff=merge_cutoff),
        expand("{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.count_ratio.csv", output_base=output_base, region=regions, cutoff=merge_cutoff),


# ------------------------------------------------------------------------------------------
# 修改后的demultiplex规则
# ------------------------------------------------------------------------------------------
rule demultiplex:
    input:
        r1 = f"{input_dir}/{{case}}_R1.fastq.gz",
        r2 = f"{input_dir}/{{case}}_R2.fastq.gz",
        primer = primer_info
    output:
        r1 = expand("{{output_base}}/demultiplex/{{case}}/{region}_R1.fastq.gz", region=regions),
        r2 = expand("{{output_base}}/demultiplex/{{case}}/{region}_R2.fastq.gz", region=regions),
    params:
        output_dir = "{output_base}/demultiplex/{case}",
        mismatch = mismatch,
        max_reads = max_reads
    log:
        "{output_base}/demultiplex/{case}.log"
    threads: 4
    shell:
        """
        python ./program/demultiplex.py \
            {input.primer} \
            {input.r1} \
            {input.r2} \
            {params.output_dir} \
            --mismatch={params.mismatch} \
            --max_reads={params.max_reads} \
            2> {log}
        """

# ------------------------------------------------------------------------------------------
# FastQC 预处理规则
# ------------------------------------------------------------------------------------------
rule fastqc:
    input:
        r1 = "{output_base}/demultiplex/{case}/{region}_R1.fastq.gz",
        r2 = "{output_base}/demultiplex/{case}/{region}_R2.fastq.gz"
    output: 
        html1 = "{output_base}/qc/fastqc/{case}/{region}_R1_fastqc.html",
        zip1 = "{output_base}/qc/fastqc/{case}/{region}_R1_fastqc.zip",
        html2 = "{output_base}/qc/fastqc/{case}/{region}_R2_fastqc.html",
        zip2 = "{output_base}/qc/fastqc/{case}/{region}_R2_fastqc.zip"
    params:
        fastqc_dir = "{output_base}/qc/fastqc/{case}"
    log:
        "{output_base}/qc/fastqc/{case}/{region}_fastqc.log"
    threads: 4
    shell:
        """
        fastqc -o {params.fastqc_dir} -t {threads} {input.r1} {input.r2} 2> {log}
        """

# ------------------------------------------------------------------------------------------
# MultiQC 规则
# ------------------------------------------------------------------------------------------
rule multiqc:
    input: 
        expand("{output_base}/qc/fastqc/{case}/{region}_R1_fastqc.zip", output_base=output_base, case=cases, region=regions),
        expand("{output_base}/qc/fastqc/{case}/{region}_R2_fastqc.zip", output_base=output_base, case=cases, region=regions)
    output:
        "{output_base}/qc/multiqc_report.html"
    params:
        fastqc_dir = "{output_base}/qc/fastqc",
        multiqc_dir = "{output_base}/qc"
    log:
        "{output_base}/qc/multiqc_report.log"
    shell: 
        """
        multiqc {params.fastqc_dir} -o {params.multiqc_dir} --filename multiqc_report.html --no-data-dir 2> {log}
        """
# ------------------------------------------------------------------------------------------
# 合并fastq文件
# ------------------------------------------------------------------------------------------
rule merge_fastq:
    input:
        r1 = "{output_base}/demultiplex/{case}/{region}_R1.fastq.gz",
        r2 = "{output_base}/demultiplex/{case}/{region}_R2.fastq.gz",
        primer = primer_info
    output:
        r1 = "{output_base}/merged/{case}/{region}_cutoff_{cutoff}_R1.fastq.gz",
        r2 = "{output_base}/merged/{case}/{region}_cutoff_{cutoff}_R2.fastq.gz",
        report = "{output_base}/merged/{case}/{region}_cutoff_{cutoff}.report"
    log:
        "{output_base}/merged/{case}/{region}_cutoff_{cutoff}.log"
    params:
        region_name = "{region}",
        cutoff = "{cutoff}"
    shell:
        """
        python ./program/merge.py \
            --input_r1 {input.r1} \
            --input_r2 {input.r2} \
            --primer_info {input.primer} \
            --region_name {params.region_name} \
            --output_r1 {output.r1} \
            --output_r2 {output.r2} \
            --report_file {output.report} \
            --cutoff {params.cutoff} > {log} 2>&1
        """
# ------------------------------------------------------------------------------------------
# 整理文件顺序
# ------------------------------------------------------------------------------------------
rule group_files:
    input:
        r1 = "{output_base}/merged/{case}/{region}_cutoff_{cutoff}_R1.fastq.gz",
        r2 = "{output_base}/merged/{case}/{region}_cutoff_{cutoff}_R2.fastq.gz",
    output:
        r1 = "{output_base}/merged_groupped_by_region/{region}/{region}_{case}_cutoff_{cutoff}_R1.fastq.gz",
        r2 = "{output_base}/merged_groupped_by_region/{region}/{region}_{case}_cutoff_{cutoff}_R2.fastq.gz"
    shell:
        """
        path_r1=`realpath {input.r1}`
        cp $path_r1 {output.r1}
        path_r2=`realpath {input.r2}`
        cp $path_r2 {output.r2}
        """
# ------------------------------------------------------------------------------------------
# 参考基因组索引文件
# ------------------------------------------------------------------------------------------
rule get_bed:
    input:
        primer_table = primer_info
    output:
        region_bed = temp('{output_base}/refs/{region}.region.bed'),
        focus_bed = temp('{output_base}/refs/{region}.focus.bed')
    shell:
        """
        python ./program/get_bed.py \
            --primer_table {input.primer_table} \
            --region_name {wildcards.region} \
            --output_region_bed {output.region_bed} \
            --output_focus_bed {output.focus_bed}
        """
rule get_fasta:
    input:
        region_bed = rules.get_bed.output.region_bed,
        focus_bed = rules.get_bed.output.focus_bed
    output:
        region_fa = temp('{output_base}/refs/{region}.ref.tmp.fa'),
        focus_fa = temp('{output_base}/refs/{region}.focus.tmp.fa')
    params:
        genome_hg38 = genome_hg38
    shell:
        """
        bedtools getfasta -nameOnly -s -bed {input.region_bed} -fi {params.genome_hg38} > {output.region_fa}
        bedtools getfasta -nameOnly -s -bed {input[1]} -fi {params.genome_hg38} > {output.focus_fa}
        """
rule rm_charater_fa:
    input:
        region_fa = rules.get_fasta.output.region_fa,
        focus_fa = rules.get_fasta.output.focus_fa
    output:
        region_fa = temp('{output_base}/refs/{region}.ref.fa'),
        focus_fa = temp('{output_base}/refs/{region}.focus.fa')
    params:
        awk = """'FNR==1{print substr($1, 1, length($1)-3);} FNR==2{print;}'"""
    shell:
        """
        awk {params.awk} {input.region_fa} > {output.region_fa}
        awk {params.awk} {input.focus_fa} > {output.focus_fa}
        """
rule letter2LETTER:
    input:
        region_fa = rules.rm_charater_fa.output.region_fa,
        focus_fa = rules.rm_charater_fa.output.focus_fa
    output:
        region_fa_upper = '{output_base}/refs/{region}.ref.upper.fa',
        focus_fa_upper = temp('{output_base}/refs/{region}.focus.upper.fa')
    params:
        awk = """'FNR==1{print $0;} FNR==2{print toupper($0);}'"""
    shell:
        """
        cat {input.region_fa} | awk {params.awk} > {output.region_fa_upper}
        cat {input.focus_fa} | awk {params.awk} > {output.focus_fa_upper}
        """
rule get_focus_seq:
    input:
        focus_fa_upper = rules.letter2LETTER.output.focus_fa_upper
    output:
        focus_seq = '{output_base}/refs/{region}.focus.upper.fa.seq'
    shell:
        """
        sed '1d' {input} > {output}
        """
rule faidx:
    input:
        '{output_base}/refs/{region}.ref.upper.fa'
    output:
        '{output_base}/refs/{region}.ref.upper.fa.fai'
    shell:
        """
        samtools faidx {input}
        """
rule bwa_idex:
    input:
        '{output_base}/refs/{region}.ref.upper.fa',
        '{output_base}/refs/{region}.ref.upper.fa.fai',
    output:
        '{output_base}/refs/{region}.ref.upper.fa.amb',
        '{output_base}/refs/{region}.ref.upper.fa.ann',
        '{output_base}/refs/{region}.ref.upper.fa.bwt',
        '{output_base}/refs/{region}.ref.upper.fa.pac',
        '{output_base}/refs/{region}.ref.upper.fa.sa'
    shell:
        """
        bwa index {input[0]}
        """
# ------------------------------------------------------------------------------------------
# 运行cutadapt
# ------------------------------------------------------------------------------------------
rule cutadapt:
    input:
        r1 = rules.group_files.output.r1,
        r2 = rules.group_files.output.r2,
        # "{output_base}/merged_groupped_by_region/{region}/{region}_{case}_cutoff_{cutoff}_R1.fastq.gz",
    output:
        r1 = temp("{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_cutadapt_R1.fastq.gz"),
        r2 = temp("{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_cutadapt_R2.fastq.gz"),
    log:
        "{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_cutadapt.log"
    params:
        cutadapt_r1 = cutadapt_r1,
        cutadapt_r2 = cutadapt_r2,
    threads: 4
    shell:
        """
        cutadapt -j {threads} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 100 -a {params.cutadapt_r1} \
        -A {params.cutadapt_r2} \
        -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1
        """
rule bwa_mapping:
    input:
        r1 = rules.cutadapt.output.r1,
        r2 = rules.cutadapt.output.r2,
        ref = rules.letter2LETTER.output.region_fa_upper
    output:
        sam = temp("{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_bwa.sam")
    log:
        "{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_bwa.log"
    threads: 4
    shell:
        """
        bwa mem {input.ref} {input.r1} {input.r2} -t {threads} -M > {output.sam} 2>{log}
        """
        # #         {BOWTIE} -x {input[2]} -1 {input[0]} -2 {input[1]} -p {THREADS} -S {output} 2>{log} not good!
rule sam_to_bam:
    input:
        sam = rules.bwa_mapping.output.sam,
        ref = rules.letter2LETTER.output.region_fa_upper
    output:
        bam = temp("{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_bwa.bam")
    shell:
        """
        samtools sort -n {input.sam} | \
        samtools view -h -f 1 -F 268 | \
        bioat bam remove_clip \
            --output_fmt BAM \
            --max_clip 6 \
            --output {output.bam}
        # -c 6 是因为使用 MGI 时会稳定有 4~6 个 softclip 在左边
        """
rule samtools_sort_by_position:
    input:
        bam = rules.sam_to_bam.output.bam,
    output:
        bam_sort = temp("{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_bwa_sort.bam")
    threads: 4
    shell:
        "samtools sort -O BAM -o {output.bam_sort} -T {output.bam_sort}.temp -@ {threads} {input.bam}"
rule samtools_index:
    input:
        bam = rules.samtools_sort_by_position.output.bam_sort
    output:
        bai = temp("{output_base}/mapping/{region}/{region}_{case}_cutoff_{cutoff}_bwa_sort.bam.bai")
    threads: 4
    shell:
        "samtools index -@ {threads} {input.bam} {output.bai}"
rule spike_in_mpileup:
    input:
        bam = rules.samtools_sort_by_position.output.bam_sort,
        bai = rules.samtools_index.output.bai,
        ref = rules.letter2LETTER.output.region_fa_upper
    output:
        mpileup = temp("{output_base}/mapping/{region}_{case}_cutoff_{cutoff}.mpileup")
    shell:
        "samtools mpileup {input.bam} --reference {input.ref} --max-depth 10000000 -q 20 -Q 20 > {output.mpileup} " 


rule parse_mpileup:
    input: 
        mpileup = rules.spike_in_mpileup.output.mpileup
    output:
        table = "{output_base}/plots/single/{region}.{case}.cutoff_{cutoff}.base_mut.tsv"
    params:
        temp_dir = "{output_base}/plots/single/{region}.{case}.cutoff_{cutoff}.base_mut.tmp"
    shell:
        "bioat bam mpileup2table {input} -o {output} --mutation_number_threshold 0 --temp_dir {params.temp_dir}"


rule bmat_plot:
    input: 
        table = rules.parse_mpileup.output.table,
        focus_seq = rules.get_focus_seq.output.focus_seq,
    output:
        pdf = "{output_base}/plots/single/{region}.{case}.cutoff_{cutoff}.base_mut.plot.pdf"
    params:
        # target_seq = lambda wildcards, input: defult_sgRNA_dict_for_plot[input[0].split("/")[4].split("-")[0]],
        target_seq = lambda wildcards, input: open(input.focus_seq, 'r').readline().strip()
    
    shell:
        """
        if [[ `cat {input.table} |wc -l` -eq 1 ]]; then 
        echo "table is empty" 
        echo "will touch a empty file" 
        touch {output.pdf}
        echo "empty! decrease the cutoff and try again" > {output.pdf}.empty.log
        else
        echo "table is ok!"
        echo "start to plot"
        bioat target_seq region_heatmap --input_table {input.table} --output_fig {output.pdf} --region_extend_length 50 --target_seq {params.target_seq} --log_level error
        fi
        """
rule bmat_multiplot:
    input: 
        table = expand("{{output_base}}/plots/single/{{region}}.{case}.cutoff_{{cutoff}}.base_mut.tsv", case=cases),
        focus_seq = rules.get_focus_seq.output.focus_seq,
        ref_seq = rules.letter2LETTER.output.region_fa_upper
    output:
        heatmap_pdf = "{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.heatmap.pdf",
        heatmap_csv = "{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.heatmap.csv",
        count_ratio_pdf = "{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.count_ratio.pdf",
        count_ratio_csv = "{output_base}/plots/multi/{region}.multiplot.cutoff_{cutoff}.base_mut.count_ratio.csv",
    params:
        input_tables = lambda wildcards, input: ",".join(input.table),
        labels = lambda wildcards: ",".join(cases),
        mut_direction = mut_direction,
        reference_seq = lambda wildcards, input: list(open(input.ref_seq, 'r').readlines())[1].strip(),
        # target_seq = lambda wildcards, input: defult_sgRNA_dict_for_plot[input[0].split("/")[4].split("-")[0]],
        target_seq = lambda wildcards, input: open(input.focus_seq, 'r').readline().strip()

    shell:
        """
        bioat target_seq region_heatmap_compare \
            --input_tables {params.input_tables} \
            --labels {params.labels} \
            --to_base A,G,C,T,Ins,Del \
            --count_ratio all \
            --heatmap_mut_direction {mut_direction} \
            --output_fig_heatmap '{output.heatmap_pdf}' \
            --output_table_heatmap '{output.heatmap_csv}' \
            --output_fig_count_ratio '{output.count_ratio_pdf}' \
            --output_table_count_ratio '{output.count_ratio_csv}' \
            --region_extend_length 0 \
            --reference_seq {params.reference_seq} \
            --target_seq {params.target_seq} \
        """
