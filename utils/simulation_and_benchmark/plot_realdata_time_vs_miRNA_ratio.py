#! /usr/bin/env python3
"""
plot the time vs miRNA ratio box plot
"""

import os, sys, re
import json
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

fn = '/data/cqs/chenh19/findadapt/findadapt_bench/2023_realdata/findadapt_res.time_mem_summary.txt'

def get_miRNA_ratio(lb):
    pw = '/data/cqs/chenh19/findadapt/findadapt_bench/2023_realdata/findadapt_res'
    fn_log = f'{pw}/{lb}.log'
    # total matched reads for inferring = 601, n_processed = 90000, matched reads density = 149.8
    
    with open(fn_log) as f:
        ratio = None
        for i in f:
            p = r'.*total matched reads for inferring = (\d+), n_processed = (\d+)'
            m = re.match(p, i)
            if m:
                matched, total_parsed = m.groups()
                total_parsed = int(total_parsed)
                matched = int(matched)
                ratio = matched / total_parsed
    
    if ratio is None:
        print(f'invalid miRNA ratio value found: {fn_log}')

    return ratio
    

fn = '/data/cqs/chenh19/findadapt/findadapt_bench/fq_size_vs_memory/findadapt_res.time_mem_summary.txt'

with open(fn) as f:
    header = f.readline().strip().split('\t')
    idx_fq = header.index('fq')
    idx_time = header.index('wall_time')
    idx_lb = header.index('lb')
    
    ratio_l = []
    ratio_d = {}
    time_l = []
    for i in f:
        line = i.strip().split('\t')
        time = line[idx_time]
        time = float(time)
        fq = line[idx_fq]

        lb = line[idx_lb]
        ratio_value = get_miRNA_ratio(lb)
        
        if not ratio_value:
            continue
        
        # so it is right open, (right border value is not included in the class)
        if ratio_value < 0.005:
            ratio_class = '0-0.5%'
        elif ratio_value < 0.01:
            ratio_class = '0.5-1%'
        elif ratio_value < 0.05:
            ratio_class = '1-5%'
        elif ratio_value < 0.1:
            ratio_class = '5-10%'
        elif ratio_value >= 0.1:
            ratio_class = '>10%'

        ratio_d[lb] = ratio_value
        ratio_l.append(ratio_class)
        time_l.append(time)

    with open(f'findadapt_fq_files_miRNA_ratio.json', 'w') as f:
        json.dump(ratio_d, f)
    

# Wall Time
order = ['0-0.5%', '0.5-1%', '1-5%', '5-10%', '>10%']
flierprops = dict(marker='o', markersize=1, linestyle='none')
fig = plt.figure(figsize=(4.5, 3),dpi=300)
a = sns.boxplot(x=ratio_l,  y=time_l, linewidth=0.4, flierprops=flierprops, order=order)
a.set_xlabel('')
a.set_ylabel('Time (s)', fontsize=12)
a.set_yscale('log')
a.set_yticks([0.1, 0.5, 1, 5, 10, 50, 100])
a.get_yaxis().set_major_formatter(ticker.ScalarFormatter())

plt.xticks(fontsize=7, rotation=30, ha='right')
plt.yticks(fontsize=10)
fig.savefig(f"findadapt_time_vs_realdata_miRNA_ratio.png", bbox_inches='tight')
