#! /usr/bin/env python3
"""
summarize the realtime runtime and smallRNA ratio
"""

import os, sys, re
from glob import glob
import pandas as pd
if __name__ == "__main__":
    pw = '/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/realdata'
    fls = glob(f'{pw}/bench_res/*.txt')
    

    res = []
    
    invalid_res = []
    for fn in fls:
        with open(fn) as f:
            lb = os.path.basename(fn).replace('findadapt.', '').replace('.txt', '')
            gse = lb.split('.', 1)[0]
            ires = [lb, gse, None, None]
            for i in f:
                try:
                    k, v = i.strip().split(':')
                    v = float(v)
                except:
                    continue
                
                if k == 'smallrna ratio':
                    if ires[2] is None or ires[2] > v:
                        ires[2] = v
                if k == 'runtime':
                    ires[3] = v
            if ires[2] is None or ires[3] is None:
                f.seek(0)
                content = f.read()
                if not 'fail to get fq adapter' in content:
                    invalid_res.append(re.match(r'.*(SR[XR]\d+)', fn).group(1))
                else:
                    print(f'fail to get adapter: {lb}')
                continue
            res.append(ires)

    if len(invalid_res) > 0:
        fn_script = f'{pw}/findadapt_cmds.sh'
        script_list = {}
        with open(fn_script) as f:
            for i in f:
                srx = re.match(r'.*(SR[XR]\d+)', i).group(1)
                script_list[srx] = i
        
        with open(f'{pw}/rerun.script.txt', 'w') as o:
            for srx in invalid_res:
                o.write(script_list[srx])


    import json
    with open(f'{pw}/realdata.benchmark.all.json', 'w') as o:
        
        json.dump(res, o)

    import json
    with open(f'{pw}/realdata.benchmark.all.json') as f:
        res = json.load(f)

    res_new = []
    time_d = {}
    for idx, i in enumerate(res):
        gse = i[1]
        time_d.setdefault(gse, []).append([idx, i[3]])  # idx and run time
    
    for v in time_d.values():
        v = sorted(v, key=lambda _: _[1])
        idx_list = [_[0] for _ in v[:2]]
        tmp = [res[_] for _ in idx_list]
        res_new.extend(tmp)
    
    def get_group(ratio):
        if ratio < 0.005:
            return '0-0.5%'

        if ratio < 0.01:
            return '0.5-1%'
        
        if ratio < 0.05:
            return '1-5%'
        
        if ratio < 0.1:
            return '5-10%'
        
        return '>=10%'
        


    data = pd.DataFrame(res_new, columns=['lb', 'gse', 'ratio', 'time'])
    data.sort_values('time', inplace=True)
    
    data['ratio_group'] = data['ratio'].map(get_group)
    
    import seaborn as sns
    import matplotlib.pyplot as plt

    sample_counts = data['ratio_group'].value_counts().sort_index()
    order = sample_counts.index


    fig = plt.figure(figsize=(6, 4),dpi=300)
    
    ax = sns.violinplot(x='ratio_group', y='time', data=data, inner=None, linewidth=0.5, order=order, color='#447eae')

    # Add stripplot to show all the dots
    sns.stripplot(x='ratio_group', y='time', data=data, color='black', size=2, jitter=True, order=order)
    
    # Add sample number to each violin plot
    for i, count in enumerate(sample_counts):
        ax.text(i, data.loc[data.ratio_group == order[i], 'time'].max(), f'n = {count}', ha='center', va='bottom', fontsize=9)

    ax.set_xlabel('')
    ax.set_ylabel('Time (s)', fontsize=12)
    
    maxtime = data['time'].max() + 1
    ax.set_ylim(0, maxtime)
    
    sns.despine()
    
    fig.savefig(f"{pw}/realdata.runtime.png", bbox_inches='tight')
