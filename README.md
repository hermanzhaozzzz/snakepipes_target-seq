# Environment:
```shell
# 使用仓库中的conda环境配置文件
conda env create -f conda_env.yml
conda activate snakepipes_target-seq
python -m ipykernel install --user --name snakepipes_target-seq # for notebook usage
```


# step.01 step.01.find_genome_region_info_for_primer_table.ipynb

该步骤核心就是完善这个表`primer_table/primer.txt` 

核心columns是： `region_index    left_primer    right_primer`

- 该文件既为文库中所涉及的靶向深度测序位置（primer amplicon region）
  - target seq region的start和end位点可以从`step.01.find_genome_region_info_for_primer_table`中获取，也可以从UCSC blat工具中获取
  - 记得 check 原始测序文件中是否有这个primer的片段，有就ok
    - (`zcat foo_R1/R2.fastq.gz | grep AAGTC... | less`)
    - 保证fwd primer是4bp random barcode，而rev primer对应6bp random barcode
  - strand信息，从Direct-seq/Detect-seq/其他测序IGV中得到：
    - CBE pattern如果是CT就是+，GA就是-
    - ABE pattern如果是AG就是+，TC就是-
    - 其他情况自行判断
      -  不确定时可以先全填`+`, 往下走流程，然后target plot出图后再修改`primer.txt` strand信息，然后重跑`step.06.form_reference.py`和`step.07.mapping_to_local_and_plot_PE.py`即可(重跑之前删除相关output: `sh clean-output_for_step6and7.sh`)



# step.02.target_seq.smk.py
- 对 config.yaml进行必要设置
- 然后运行整个项目

```
snakemake -p -j 50 -s step.02.demultiplex_fastq.smk.py --rerun-incomplete --cores 8 -n
```
# step.03.complicated_analysis.ipynb
- 顾名思义,后期处理的一些示例代码