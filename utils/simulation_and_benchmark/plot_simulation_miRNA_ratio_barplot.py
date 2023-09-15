#! /usr/bin/env python3
"""
bar plot for the simulation data for different miRNA ratio, memory and time usage
"""
import seaborn as sns
import matplotlib.pyplot as plt

fn = '/data/cqs/chenh19/findadapt/findadapt_bench/simu_mirna_ratio/findadapt_res.time_mem_summary.txt'

with open(fn) as f:
    header = f.readline().strip().split('\t')
    idx_fq = header.index('fq')
    idx_lb = header.index('lb')
    idx_time = header.index('wall_time')
    idx_mem = header.index('max_mem')
    
    
    time_l = []
    mem_l = []
    lb_l = []
    
    for i in f:
        line = i.strip().split('\t')
        lb = line[idx_lb]
        time = line[idx_time]
        mem = line[idx_mem]
        lb = lb.rsplit('_', 1)[-1]

        # if lb == '0.0001':
        #     continue
        
        time_l.append(float(time))
        mem_l.append(float(mem))
        
        # benchmark_test.miRNA_ratio_0.01
        lb_l.append(lb)


pwout = '/data/cqs/chenh19/findadapt/findadapt_bench/simu_mirna_ratio'
# plot
fig = plt.figure(figsize=(4.5, 3),dpi=300)
a = sns.barplot(x=lb_l,  y=time_l, color='grey')
a.set_xlabel('miRNA ratio', fontsize=12)
a.set_ylabel('Time (s)', fontsize=12)
# plt.xticks(fontsize=7, rotation=30, ha='right')
plt.yticks(fontsize=10)
fig.savefig(f"{pwout}/miRNA_ratio_vs_time.png", bbox_inches='tight')


fig = plt.figure(figsize=(4.5, 3),dpi=300)
a = sns.barplot(x=lb_l,  y=mem_l, color='grey')
a.set_xlabel('miRNA ratio', fontsize=12)
a.set_ylabel('Max Memory (MB)', fontsize=12)
# plt.xticks(fontsize=7, rotation=30, ha='right')
plt.yticks(fontsize=10)
fig.savefig(f"{pwout}/miRNA_ratio_vs_mem.png", bbox_inches='tight')


