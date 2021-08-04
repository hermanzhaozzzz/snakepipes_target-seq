#!/bin/sh
# properties = {"type": "single", "rule": "MergeFastq", "local": false, "input": ["../TargetSeq-CBE-HEK4-rep2"], "output": ["../TargetSeq-CBE-HEK4-rep2/cutoff_5/CBE-HEK4-rep2_cutoff_5.report"], "wildcards": {"lib": "CBE-HEK4-rep2", "cutoff": "5"}, "params": {"out_dir": "../TargetSeq-CBE-HEK4-rep2/cutoff_5", "cutoff": "5"}, "log": ["../TargetSeq-CBE-HEK4-rep2/cutoff_5/CBE-HEK4-rep2_cutoff_5.log"], "threads": 1, "resources": {}, "jobid": 23, "cluster": {}}
 cd /lustre1/chengqiyi_pkuhpc/zhaohn/3.project/2021_ABE_topic/20210727_TargetSeq_batch1_ABE_ABE8e_ACBE_CBE-HEK4/snakepipes_target-seq && \
PATH='/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/miniconda3/bin':$PATH /lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/miniconda3/bin/python3.8 \
-m snakemake ../TargetSeq-CBE-HEK4-rep2/cutoff_5/CBE-HEK4-rep2_cutoff_5.report --snakefile /lustre1/chengqiyi_pkuhpc/zhaohn/3.project/2021_ABE_topic/20210727_TargetSeq_batch1_ABE_ABE8e_ACBE_CBE-HEK4/snakepipes_target-seq/step.04.merge_reads_for_regions.py \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /lustre1/chengqiyi_pkuhpc/zhaohn/3.project/2021_ABE_topic/20210727_TargetSeq_batch1_ABE_ABE8e_ACBE_CBE-HEK4/snakepipes_target-seq/.snakemake/tmp.zo4mmf_y ../TargetSeq-CBE-HEK4-rep2 --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules MergeFastq --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /lustre1/chengqiyi_pkuhpc/zhaohn/3.project/2021_ABE_topic/20210727_TargetSeq_batch1_ABE_ABE8e_ACBE_CBE-HEK4/snakepipes_target-seq/.snakemake/tmp.zo4mmf_y/23.jobfinished || (touch /lustre1/chengqiyi_pkuhpc/zhaohn/3.project/2021_ABE_topic/20210727_TargetSeq_batch1_ABE_ABE8e_ACBE_CBE-HEK4/snakepipes_target-seq/.snakemake/tmp.zo4mmf_y/23.jobfailed; exit 1)

