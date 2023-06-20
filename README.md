# Environment:
```shell
# 使用仓库中的conda环境配置文件
conda env create -f conda_env.yml
# 由于脚本升级暂未完成，同时依赖一个python2.7环境，如手动创建环境
conda create -n snakepipes_py27 \
    python=2.7.15 \
    biopython=1.72 \
    pandas=0.24.2 \
    numpy=1.16.5 \
    matplotlib \
    bedtools=2.29.2 \
    bwa \
    samtools
```

# 先进行整体文库质控
整体质控先跑一次，直接使用之前写的snakepipes
https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc
```shell
git clone git@github.com:hermanzhaozzzz/snakepipes_fastqc-multiqc.git

# 接下来按照该仓库要求进行质控即可
```

# step.01 `step.01.find_genome_region_info_for_primer_table.ipynb`

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






# step.02 split raw FASTQ
- 接下来跑脚本按照primer info拆分不同的target位点
- Python2的`[TODO] 这里打算整合到bioat里`
- 这一步是脚本一次性扔到后台去跑

```
bash step.02.sh.01.demultiplex_fastq_pipline.sh
```
```
tail -f ../TargetSeq*/*.log
# 北极星会将日志重写到stderr
tail -f *.err
```



# step.03 sampling raw reads
```
# 直接取一定量的测序结果继续分析
# 也可把测序数据全用上，但没必要，比较耗时和耗内存
# target deep sequencing 达到一定深度即可

# local
snakemake -pr -j 100 -s step.03.head-reads_for_each-region.py -n
# polaris (根据内存和核心数目适当调用队列、核心数目、任务数量)
snakemake --cluster "pkurun-cnlong 1 4" -pr -j 100 -s step.03.head-reads_for_each-region.py -n

```
# step.04 单位点质控

```shell
# local
snakemake -pr -j 15 -s step.04.qc_for_regions.py -n
# polaris
snakemake --cluster 'pkurun-cnnl 1 20' -pr -j 15 -s step.04.qc_for_regions.py -n
```

在`../qc_quantify-all/multiqc/multiqc_report.html`查看单位点质量，双端raw reads总量不足1.5M提示可能实验无效或者不准确，结合最终有效mapping reads数进行判断

# step.05 UMI 校正测序错误和PCR错误
- 这一步是通过barcode来merge duplicate和校正测序错误、PCR错误
- 注意选择cutoff！默认选`cutoff=3`
- check一下cutoff取值的区别【cutoff越高越严格，意思是几个reads互相校正】
```
# local 注意内存开销
snakemake -pr -j 2 step.05.merge_reads_for_regions.py -n
tail -f ../TargetSeq*/cutoff*/*.log
# polaris
snakemake --cluster 'pkurun-cnnl 1 20' -pr -j 6 step.05.merge_reads_for_regions.py -n
tail -f ../TargetSeq*/cutoff*/*.log
```
# step.06 生成扩增位点reference信息
```
# 不耗资源直接login节点或者本地运行即可
# local
snakemake -pr -j 100 -s step.06.form_reference.py -n
```

# step.07 计算单碱基高精度突变信息并绘制target plot

```
# local
snakemake -pr -j 10 -s step.07.mapping_to_local_and_plot_PE.py -n

# polaris
snakemake --cluster "pkurun-cnlong 1 4" -pr -j 100 -s step.07.mapping_to_local_and_plot_PE.py -n

```



`step.07.mapping_to_local_and_plot_SE.py`为1.5代测序所用，没用到时不需要了解

# step.08 绘制multi-treatment target plot

`详见README.ipynb`
