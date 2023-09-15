import matplotlib.pyplot as plt
import pandas as pd

fn = '/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/realdata/87_geo_study_group.txt'

reads_type = ['3p Adapter only', 'QIAseq', 'NEXTFLEX', 'Template Switching', 'Trimmed']

ct = {k: [0 for _ in range(5)] for k in ['Matched', 'Mislabeled', 'Not Specified', 'Other']}

idx_map = {v: n for n, v in enumerate(reads_type)}

with open(fn) as f:
    f.readline()
    for i in f:
        group, subgroup = i.strip().split('\t')
        if subgroup == 'other':
            subgroup = 'Other'
        idx = idx_map[group]
        
        ct[subgroup][idx] += 1
        
        
plotdata = pd.DataFrame(ct, index = reads_type)

pwout = '/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench'



# plot
ax = plotdata.plot(kind='bar', stacked = True, figsize=(4, 3), color=['green', 'orange', 'lightblue', 'grey'])

fontsize = 6
# plt.legend(loc='upper left', bbox_to_anchor=(1.01, 0.7), fontsize=6)

ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=fontsize)
# ax.set_yticklabels(ax.get_yticklabels(), fontsize=fontsize)
ax.set_ylabel('Number of studies')
legend = ax.get_legend()

# legend._legend_box = None
# legend._legend_handle_box = None
legend.borderpad = 1
legend.labelspacing = 0
legend.handlelength = 0
legend.handletextpad = 0
for text in legend.get_texts():
    text.set_fontsize(fontsize)
# plt.setp(legend.get_texts(), fontsize=fontsize) 

plt.savefig(f"{pwout}/87_geo_acuracy.png", dpi=300, bbox_inches='tight')



# handles, labels = ax.get_legend_handles_labels()
# filtered_handles_labels = [(handle, label) for handle, label in zip(handles, labels) if label != 'other']
# ax.legend([handle for handle, _ in filtered_handles_labels], [label for _, label in filtered_handles_labels])
