# 1. primer design table -> primer amplicon region
- TSS和TES可以从UCSC genome browser  -> Tools -> Blat查询比对
    ![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh5nc9yeutj30mc0avjwl.jpg)
    https://genome.ucsc.edu/cgi-bin/hgBlat?command=start
    ![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh5ozmlaurj31jg0u0av2.jpg)
    然后将./primer_table里的xlsx转成txt格式和
- 找到strand
    方向！找strand
    T for +
    A for -
# 2. split raw FASTQ
- 接下来跑脚本按照primer info拆分不同的target位点
- Python2的
- 这一步！是脚本一次性扔到后台去跑，但是完成了脚本要htop监控或查看lig看一下到底跑完没

```
bash sh.01.demultiplex_fastq_pipline.sh  
```
```
tail -f ../TargetSeq*/*.log           
```

default reads =1 000 000

# 3. merge and validation of FASTQ
- 这一步是通过barcode来merge duplicate和校正测序错误
- 注意选择cutoff！
```
snakemake --jobs 6 --cores 24 --snakefile snk.03.merge_fastq_pipline.py -p
```

```
tail -f ../TargetSeq*/cutoff*/*.log           
```
# 4. get ref fasta
```
snakemake --jobs 6 --cores 24 --snakefile snk.04.getfast.py -p
```

# 5. plot base mutation stat by APOBEC
```
snakemake --jobs 6 --cores 24 --snakefile snk.05.merge_fastq_to_plot_pipline.bwa.py -p
```


1. 环境
snakemake                 5.3.0
python                    2.7.15
biopython                 1.76  !1.72
numpy                     1.16.5
pandas                    0.24.2 
bedtools                  2.29.2

2. bug1: 327行 【.】 得分为1？？？
![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh9bi7q8vmj30s40ilq49.jpg)
怀疑是第一步拆分没有跑完的问题，跑完再check一下
3. bug2： 855行，添加了 for 循环，解决 x没有定义的问题
4. plot的字体问题
没有Arial，用了Default
5. bug 点的问题
line 312
![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh9iqqqhidj30y80mkq7c.jpg)

![](https://tva1.sinaimg.cn/large/007S8ZIlly1gh9iww7avoj30iw08mmxk.jpg)



全转大写字母
UCSC标记为准
getfasta的index问题
