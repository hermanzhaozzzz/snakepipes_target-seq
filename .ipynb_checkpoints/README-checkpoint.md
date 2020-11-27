# Env:
```
conda create -n snakepipes_target-seq-from-table-to-plot python=2.7.15 biopython=1.72 pandas=0.24.2 numpy=1.16.5 matplotlib bedtools=2.29.2 bwa
```

# 1. primer design table -> primer amplicon region
target seq region的start和end位点可以从step.01的ipynb protocol中获取，也可以从UCSC blat工具中获取

记得check原始测序文件中是否有这个primer的片段，有就ok

strand信息：
    如果知道sgRNA的话就是该sgRNA的strand方向
    CBE pattern如果是CT就是+，GA就是-
    不确定时可以先往下走流程，然后target plot出图后再修改第5、6步即可
```
zcat foo.fq.gz | grep primer-sequence
```
    
# 2. 先跑质控的pipeline
整体质控先跑一次，直接使用之前写的snakepipes
https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc
```
git@github.com:hermanzhaozzzz/snakepipes_fastqc-multiqc.git
```




# 3.split raw FASTQ，跑step2
- 接下来跑脚本按照primer info拆分不同的target位点
- Python2的
- 这一步是脚本一次性扔到后台去跑，但是完成了脚本要htop监控或查看log看一下到底跑完没，比较慢，1~3h

这个alloc节点之后跑，sh01-1单独占用node，sh01-2~6单独占用node sh01-7~8单独占用node，共三个node
```
bash step.02.sh.01.demultiplex_fastq_pipline.sh
```
```
tail -f ../TargetSeq*/*.log           
```



# 4. 拆分后的质控，跑step3
```
snakemake -pr -j 48 -s step.03.snk.quantify_for_every_region.py -n
```
# 3. merge and validation of FASTQ， 跑step4
- 这一步是通过barcode来merge duplicate和校正测序错误
- 注意选择cutoff！
- check一下cutoff取值的区别【cutoff越高越严格，意思是几个reads互相校正】
## 注意和建议！！！！！
用snakemake --profile slurm参数来运行任务，不容易导致崩溃】，因为这样子一个node只跑一个，我测试过，如果salloc node跑任务，一个node最多跑三个merge的并行任务，要不然会因为内存溢出而报错
```
snakemake --profile slurm -pr -j 6 -s step.04.snk.merge_fastq_pipline.py -n

tail -f ../TargetSeq*/cutoff*/*.log           
```
# 4. get ref fasta，跑step5
先使用flag_for_snk04_and_snk05.ipynb这个笔记本文件得到需要的SAMPLE list，并更改snk04的文件内容
再运行
直接提交在登录节点即可，耗资源很少，而且slurm队列IO很慢反而降低效率
注意两点
1. 记得更改txt路径信息！
```
snakemake -pr -j 48 -s step.05.snk.getfasta.py -n
```

# 5. plot base mutation stat， step6

```
snakemake -pr -j 48 -s step.06.snk.merge_fastq_to_plot_pipline.bwa.py -n
```


# 6. multiplot，step7，详见README.ipynb
