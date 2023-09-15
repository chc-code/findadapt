#! /usr/bin/env python3
"""
build the bar plot for simulation accuracy
"""

import os, sys, re
import json
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

pwout = '/nobackup/h_vangard_1/chenh19/findadapt_bench'

if __name__ == "__main__":
    fn = f'{pwout}/bench_res/benchmark_final_res.json'
    
    with open(fn) as f:
        data = json.load(f)
        # k1 = simu_type
        # k2 = tool name
        # v = list of 20 element,
        # ele = 0/1,  1=result is correct, 0 = result = wrong
        
        res = []
        for simu_type, v1 in data.items():
            for toolname, v2 in v1.items():
                score = sum(v2)
                if len(v2) != 20:
                    print(f'error, simu times != 20: {simu_type}:  {toolname}')
                    continue
                ratio = score / len(v2)
                res.append([simu_type, toolname, ratio])
        
                
        
    df = pd.DataFrame(res, columns=['simu_type', 'tool', 'accuracy'])
    
    hue_order =   [
            'findadapt',
            'dnapi',
            'earrings_hg38_sensitive',
            'earrings_hg38',
            'atropos_kmer',
            'atropos_heuristic',
            'fastp',
    ]
    
    hue_name_formal = [
            'FindAdapt',
            'DNApi',
            'EARRINGS sensitive',
            'EARRINGS',
            'Atropos khmer',
            'Atropos heuristic',
            'fastp',
    ]
    
    # order x axis
    simu_type_order = {
        'add_5mer': [1, "Adapter + 5' random-mer"],
        'add_3mer': [2, "Adapter + 3' random-mer"],
        'add_both_side': [3, "Adapter + 5' and 3' random-mer"],

        'polya_add_5mer': [4, "Poly(A) + 5' random-mer"],
        'polya_add_3mer':[5, "Poly(A) + 3' random-mer"],
        'polya_add_both_side': [6, "Poly(A) + 5' and 3' random-mer"],
        
        }

    df_plot = df.loc[df['simu_type'].isin(simu_type_order)]

    df_plot['order'] = df_plot['simu_type'].apply(lambda _: simu_type_order[_][0] - 1)
    df_plot['simu_type_new'] = df_plot['simu_type'].apply(lambda _: simu_type_order[_][1])
    df_plot = df_plot.sort_values(['order']).reset_index(drop=True)

    # plot
    fig = plt.figure(figsize=(4.5, 3),dpi=300)
    a = sns.barplot(x='simu_type_new',  y='accuracy', hue='tool', hue_order=hue_order, data=df_plot)
    # a.axes.set_title(f'Simulation Datasets Accuracy', fontsize=12)
    a.set_xlabel('')
    a.set_ylabel('Accuracy', fontsize=12)
    plt.xticks(fontsize=7, rotation=30, ha='right')
    plt.yticks(fontsize=10)

    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 0.7), fontsize=6)
    legend = plt.gca().get_legend()
    for t, l in zip(legend.texts, hue_name_formal):
        t.set_text(l)

    fig.savefig(f"{pwout}/bench_res/reads_strunc_simulation_with_randomer.png", bbox_inches='tight')

    # do the adapter only one separatedly
    # 'adapter_only': [0, "Adapter Only"],
    #    'polya_only': [4, 'Poly(A)'],
    #    'no_adapt': [7, 'No adapter'],

    simu_type_order = {
        'adapter_only': [0, "Adapter Only"],
        'polya_only': [1, 'Poly(A)'],
        'no_adapt': [2, 'No adapter'],
        }
        
    df_plot = df.loc[df['simu_type'].isin(simu_type_order)]
    df_plot['order'] = df_plot['simu_type'].apply(lambda _: simu_type_order[_][0] - 1)
    df_plot['simu_type_new'] = df_plot['simu_type'].apply(lambda _: simu_type_order[_][1])
    df_plot = df_plot.sort_values(['order']).reset_index(drop=True)

    # plot
    fig = plt.figure(figsize=(4.5, 3),dpi=300)
    a = sns.barplot(x='simu_type_new',  y='accuracy', hue='tool', hue_order=hue_order, data=df_plot)
    # a.axes.set_title(f'Simulation Datasets Accuracy', fontsize=12)
    a.set_xlabel('')
    a.set_ylabel('Accuracy', fontsize=12)
    plt.xticks(fontsize=7, rotation=30, ha='right')
    plt.yticks(fontsize=10)

    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 0.7), fontsize=6)
    legend = plt.gca().get_legend()
    for t, l in zip(legend.texts, hue_name_formal):
        t.set_text(l)

    fig.savefig(f"{pwout}/bench_res/reads_strunc_simulation.adapter_only.png", bbox_inches='tight')
