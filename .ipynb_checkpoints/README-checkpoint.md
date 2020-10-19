# 1. primer design table -> primer amplicon region
- TSS和TES可以从UCSC genome browser  -> Tools -> Blat查询比对
    ![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh5nc9yeutj30mc0avjwl.jpg)
    https://genome.ucsc.edu/cgi-bin/hgBlat?command=start
    ![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh5ozmlaurj31jg0u0av2.jpg)
    然后将./primer_table里的xlsx转成txt格式和
- check primer 
```
zcat foo.fq.gz | grep primer-sequence
```
- 找到strand
    方向！从IGV里找strand
    T for +
    A for -
    
    

    
    
    
    
# 2. 先跑质控的pipeline
直接使用之前写的snakepipes
https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc
```
git@github.com:hermanzhaozzzz/snakepipes_fastqc-multiqc.git
```




# 3. split raw FASTQ
- 接下来跑脚本按照primer info拆分不同的target位点
- Python2的
- 这一步是脚本一次性扔到后台去跑，但是完成了脚本要htop监控或查看log看一下到底跑完没，比较慢，1~3h

这个alloc节点之后跑，sh01-1单独占用node，sh01-2~6单独占用node sh01-7~8单独占用node，共三个node
```
bash sh.01.demultiplex_fastq_pipline.sh  
```
```
tail -f ../TargetSeq*/*.log           
```

# 3. merge and validation of FASTQ
- 这一步是通过barcode来merge duplicate和校正测序错误
- 注意选择cutoff！
- check一下cutoff取值的区别【还未测试】
## 注意和建议！！！！！
用snakemake --profile slurm参数来运行任务，不容易导致崩溃】，因为这样子一个node只跑一个，我测试过，如果salloc node跑任务，一个node最多跑三个merge的并行任务，要不然会因为内存溢出而报错
```
snakemake --profile slurm --jobs 6 --snakefile snk.03.merge_fastq_pipline.py -p
```

```
tail -f ../TargetSeq*/cutoff*/*.log           
```
# 4. get ref fasta
先使用flag_for_snk04_and_snk05.ipynb这个笔记本文件得到需要的SAMPLE list，并更改snk04的文件内容
再运行
直接提交在登录节点即可，耗资源很少，而且slurm队列IO很慢反而降低效率
注意两点
1. 记得更改txt路径信息！
2. 多跑几次，容易漏
```
snakemake --jobs 6 --cores 24 --snakefile snk.04.getfast.py -p
```

# 5. plot base mutation stat by APOBEC
同样先使用flag_for_snk04_and_snk05.ipynb更改SAMPLE list再依次运行snk05各个子文件
我对snk05做了一个优化，使得可以在rule plot的时候自动识别sgRNA的gene site，和自动选择适合的sgRNA sequence
![](https://tva1.sinaimg.cn/large/007S8ZIlly1ghw1xhg2tvj31f00u043m.jpg)
![](https://tva1.sinaimg.cn/large/007S8ZIlly1ghw1y7u9vrj31n00huwig.jpg)


这个alloc节点之后跑，可以跑满24个核心，占用内存资源少，占用较多计算资源
```
snakemake --jobs 6 --cores 24 --snakefile snk.05.merge_fastq_to_plot_pipline.bwa.py -p
```


1. 环境
```
conda create -n snakepipes_target-seq-from-table-to-plot python=2.7.15 biopython=1.72 pandas=0.24.2 numpy=1.16.5 matplotlib bedtools=2.29.2 bwa
```


2. BUGs:
- plot的字体问题
    没有Arial，用了Default
- 点的问题:
    脚本line 312
    是Biopython版本问题，已解决

![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh9iqqqhidj30y80mkq7c.jpg)

![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh9iww7avoj30iw08mmxk.jpg)




# 6. check igv
```
$ tar -zcvf igv.tgz TargetSeq*/cutoff_3/map*/*sort.bam* TargetSeq_BED_sample_lib
```
