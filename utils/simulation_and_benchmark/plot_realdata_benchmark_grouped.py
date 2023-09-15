#! /usr/bin/env python3
"""
build the bar plot for realdata accuracy
"""

import os, sys, re
import json
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

pw = '/data/cqs/chenh19/findadapt/findadapt_bench/2023_realdata'

def summarize_data():
    tools = ['findadapt', 'DNApi', 'EARRINGS', 'EARRINGS_sensitive', 'atropos_heuristic', 'atropos_kmer', 'fastp']

        # lb	fq	wall_time	max_mem	gse	kit_exp	seq_3p_exp	seq_3p	phase_3p	ct_3p	phase_5p	ct_5p	match
        # SRR5755813	/data/cqs/chenh19/findadapt/findadapt_bench/download/SRR5755813.fastq.gz	0.17	16.47265625	GSE100467	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG	TGGAATTCTCGG	0	783	0	989	match


    res = {}  # k1 = adapter/kit type, k2 = tool,  v = [lb, time, memory, match_result]
    for itool in tools:
        fn_summary = f'{pw}/{itool}_res.time_mem_summary.txt'
        with open(fn_summary) as f:
            header = f.readline().strip().split('\t')
            colmap = {k: idx for idx, k in enumerate(header)}
            for line in f:
                line = line.strip().split('\t')
                kit = line[colmap['kit_exp']]
                seq_exp = line[colmap['seq_3p_exp']]
                match = line[colmap['match']]
                match = 1 if 'match' in match else 0
                time = float(line[colmap['wall_time']])
                max_mem = float(line[colmap['max_mem']])
                lb = line[colmap['lb']]
                if 'NNNN' in seq_exp:
                    adapter_type = 'polyN'
                elif ' ' in seq_exp:
                    adapter_type = 'Other'
                elif 'Truseq' in kit:
                    adapter_type = 'TruSeq'
                elif seq_exp == 'empty':
                    adapter_type = 'Trimmed'
                elif kit == 'NA' or 'smallRNA v1.5' in kit:
                    adapter_type = 'Other'
                else:
                    adapter_type = kit
                res.setdefault(adapter_type, {}).setdefault(itool, []).append([lb, time, max_mem, match])

    l = []
    accuracy = []
    for adapter_type, v1 in res.items():
        for itool, v2 in v1.items():
            match_l = [_[-1] for _ in v2]
            tmp = sum(match_l) / len(match_l)
            accuracy.append([adapter_type, itool, tmp])
            for i in v2:
                l.append([adapter_type, itool] + i)

    
    
    df = pd.DataFrame(l, columns=['adapter_type', 'tool', 'SRX', 'wall_time', 'max_mem', 'match'])
    df.to_csv('studies_without_random_seq.combined.for_plot.csv', index=False)
    
    df_accuracy = pd.DataFrame(accuracy, columns=['adapter_type', 'tool', 'accuracy'])
    df_accuracy.to_csv('studies_without_random_seq.combined.accuracy.csv', index=False)
    
    return df, df_accuracy
    
    


if __name__ == "__main__":

        
    df, df_accuracy = summarize_data()
    
    hue_order =   ['findadapt', 'DNApi', 'EARRINGS', 'EARRINGS_sensitive', 'atropos_heuristic', 'atropos_kmer', 'fastp']
    
    hue_name_formal = [
            'FindAdapt',
            'DNApi',
            'EARRINGS',
            'EARRINGS sensitive',
            'Atropos heuristic',
            'Atropos khmer',
            'fastp',
    ]
    
    # order x axis
    adapter_type_order = {k: idx for idx, k in enumerate(['TruSeq', 'NEBNext', 'QIAseq', 'polyN', 'Trimmed', 'Other'])}

    df_plot = df_accuracy.loc[df_accuracy['adapter_type'].isin(adapter_type_order)]
    df_plot['order'] = df_plot['adapter_type'].apply(lambda _: adapter_type_order[_])
    df_plot = df_plot.sort_values(['order']).reset_index(drop=True)

    # plot
    fig = plt.figure(figsize=(4.5, 3),dpi=300)
    a = sns.barplot(x='adapter_type',  y='accuracy', hue='tool', hue_order=hue_order, data=df_plot)
    # a.axes.set_title(f'Simulation Datasets Accuracy', fontsize=12)
    a.set_xlabel('')
    a.set_ylabel('Accuracy', fontsize=12)
    plt.xticks(fontsize=7, rotation=30, ha='right')
    plt.yticks(fontsize=10)
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 0.7), fontsize=6)
    legend = plt.gca().get_legend()
    for t, l in zip(legend.texts, hue_name_formal):
        t.set_text(l)
    fig.savefig(f"{pw}/studies_without_random_seq.accuracy.png", bbox_inches='tight')

    # do the box plot for memory usage and time
    

    df_plot = df.loc[df['adapter_type'].isin(adapter_type_order)]
    df_plot['order'] = df_plot['adapter_type'].apply(lambda _: adapter_type_order[_])
    df_plot = df_plot.sort_values(['order']).reset_index(drop=True)
    
    
    # Wall Time
    flierprops = dict(marker='o', markersize=1, linestyle='none')
    fig = plt.figure(figsize=(4.5, 3),dpi=300)
    a = sns.boxplot(x='adapter_type',  y='wall_time', hue='tool', hue_order=hue_order, data=df_plot, linewidth=0.4, flierprops=flierprops)
    a.set_xlabel('')
    a.set_ylabel('Time (s)', fontsize=12)
    a.set_yscale('log')
    a.set_yticks([0.1, 0.5, 1, 5, 10, 50, 100, 500])
    a.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    
    plt.xticks(fontsize=7, rotation=30, ha='right')
    plt.yticks(fontsize=10)
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 0.7), fontsize=6)
    legend = plt.gca().get_legend()
    for t, l in zip(legend.texts, hue_name_formal):
        t.set_text(l)
    fig.savefig(f"{pw}/studies_without_random_seq.walltime.png", bbox_inches='tight')


    # Max Memroy
    flierprops = dict(marker='o', markersize=1, linestyle='none')
    fig = plt.figure(figsize=(4.5, 3),dpi=300)
    a = sns.boxplot(x='adapter_type',  y='max_mem', hue='tool', hue_order=hue_order, data=df_plot, linewidth=0.4, flierprops=flierprops)
    a.set_xlabel('')
    a.set_ylabel('Max Memory Usage (MB)', fontsize=12)
    a.set_yscale('log')
    a.set_yticks([10, 50, 100, 500, 1000, 5000])
    a.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    
    plt.xticks(fontsize=7, rotation=30, ha='right')
    plt.yticks(fontsize=10)
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 0.7), fontsize=6)
    legend = plt.gca().get_legend()
    for t, l in zip(legend.texts, hue_name_formal):
        t.set_text(l)
    fig.savefig(f"{pw}/studies_without_random_seq.max_mem.png", bbox_inches='tight')


