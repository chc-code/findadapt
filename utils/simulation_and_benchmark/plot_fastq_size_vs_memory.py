#! /usr/bin/env python3
"""
plot the memroy usage vs the fastq file size
"""
import os
import seaborn as sns
import matplotlib.pyplot as plt

fn = '/data/cqs/chenh19/findadapt/findadapt_bench/fq_size_vs_memory/findadapt_res.time_mem_summary.txt'

with open(fn) as f:
    header = f.readline().strip().split('\t')
    idx_fq = header.index('fq')
    idx_mem = header.index('max_mem')
    
    mem_l = []
    size_l = []
    
    for i in f:
        line = i.strip().split('\t')
        fn_fq = line[idx_fq]
        mem = line[idx_mem]
        
        size = os.path.getsize(fn_fq) / 1024 / 1024

        size_l.append(size)
        mem_l.append(float(mem))

pwout = '/data/cqs/chenh19/findadapt/findadapt_bench/fq_size_vs_memory'
# plot
max_mem = max(mem_l) + 1
fig = plt.figure(figsize=(4.5, 3),dpi=300)
plt.ylim(0, max_mem)
a = sns.scatterplot(x=size_l,  y=mem_l, s=20, edgecolors='none')
a.set_xlabel('file size (MB)', fontsize=12)
a.set_ylabel('Max Memory (MB)', fontsize=12)

# plt.xticks(fontsize=7, rotation=30, ha='right')
plt.yticks(fontsize=10)
fig.savefig(f"{pwout}/memory_vs_fq_size.png", bbox_inches='tight')


