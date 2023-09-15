#! /usr/bin/env python3
"""
bench mark the real datasets
"""
import os, sys, re
import json
import platform
node = platform.node()
if 'bioinfo2' in node:
    pw = '/data/findadapt'
    node = 'bioinfo2'
elif 'hanuman' in node:
    pw = '/data/cqs/chenh19/findadapt/findadapt_bench'
    node = 'hanuman'

study_list = """
GSE106303	QIAseq	AACTGTAGGCAC
GSE123473	NEBNext	AGATCGGAAGAG
GSE129255	NEBNext	AGATCGGAAGAG
GSE130316	NEBNext	AGATCGGAAGAG
GSE130654	NEBNext	AGATCGGAAGAG
GSE138983	NEBNext	AGATCGGAAGAG
GSE174224	NEBNext	AGATCGGAAGAG
GSE185273	NEBNext	AGATCGGAAGAG
GSE85761	NEBNext	AGATCGGAAGAG
GSE50429	NA	ATCTCGTATGCC
GSE149645	NA	empty
GSE163507	NA	empty
GSE106277	NA	NNNNNNNNNNNN
GSE112289	NA	empty
GSE128803	NA	empty
GSE141326	NA	empty
GSE150460	NA	empty
GSE97644	NA	empty
GSE144781	NA	empty
GSE67004	NA	empty
GSE138107	NA	empty
GSE117479	NA	empty
GSE126049	NA	empty
GSE128004	NA	empty
GSE134205	NA	empty
GSE113411	NEBNext	AGATCGGAAGAG
GSE156380	NEBNext and Truseq/TailorMix/Lexogen/CleanTag	AGATCGGAAGAG TGGAATTCTCGG
GSE45722	Truseq/TailorMix/Lexogen/CleanTag	AGATCGGAAGAG TGGAATTCTCGG
GSE117744	NEBNext	AGATCGGAAGAG
GSE114008	NA	NNNNNNNNNNNN
GSE146880	NA	GATCGGAAGAGC
GSE106238	NEBNext	AGATCGGAAGAG
GSE109356	NEBNext	AGATCGGAAGAG
GSE111803	NEBNext	AGATCGGAAGAG
GSE122488	NEBNext	AGATCGGAAGAG
GSE123972	NEBNext	AGATCGGAAGAG
GSE125905	NEBNext	AGATCGGAAGAG
GSE128359	NEBNext	AGATCGGAAGAG empty
GSE136321	NEBNext	AGATCGGAAGAG
GSE140106	NEBNext	AGATCGGAAGAG
GSE143613	NEBNext	AGATCGGAAGAG
GSE144627	NEBNext	AGATCGGAAGAG
GSE147517	NEBNext	AGATCGGAAGAG
GSE155281	NEBNext	AGATCGGAAGAG
GSE157528	NEBNext	AGATCGGAAGAG
GSE158659	NEBNext	AGATCGGAAGAG
GSE161339	NEBNext	AGATCGGAAGAG
GSE165323	NEBNext	AGATCGGAAGAG
GSE176134	NEBNext	AGATCGGAAGAG
GSE53451	NEBNext	AGATCGGAAGAG
GSE58410	NEBNext	AGATCGGAAGAG
GSE71008	NEBNext	AGATCGGAAGAG
GSE71661	NEBNext	AGATCGGAAGAG
GSE83539	NEBNext	AGATCGGAAGAG
GSE83669	NEBNext	AGATCGGAAGAG
GSE84306	NEBNext	AGATCGGAAGAG
GSE93175	NEBNext	AGATCGGAAGAG
GSE99430	NEBNext	AGATCGGAAGAG
GSE134220	QIAseq	AACTGTAGGCAC
GSE167863	QIAseq	AACTGTAGGCAC
GSE182717	QIAseq	AACTGTAGGCAC
GSE38916	smallRNA v1.5(illumina)	ATCTCGTATGCC
GSE100467	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE103831	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE106453	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE115572	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE122656	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE136997	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE142819	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE155360	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE179382	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE70432	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE74759	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE93020	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
GSE94721	Truseq/TailorMix/Lexogen/CleanTag	TGGAATTCTCGG
""".split('\n')
study_list = [_.strip().split('\t') for _ in study_list if _.strip()]


def get_fq_path(fn_list):
    fls = []
    not_found = []
    with open(fn_list) as f:
        for i in f:
            i = i.strip()
            if not i:
                continue
            _, fn = i.split('\t')
            fn_pure = os.path.basename(fn)
            fls.append(fn)
    fls = sorted(fls)
    return fls
            

def check_fastq_fls():
    ct = {}
    res = []
    res_lite = []
    
    n_fq_per_study = 5
    meta = {}
    not_found = []
    fq_list_for_testing = []
    for gse, kit, seq in study_list:
        # the seq can be multiple adapters
        
        fn_fqlist = f'{pw}/fqlist/{gse}.fq_list.txt'
        if not os.path.exists(fn_fqlist):
            ct[gse] = {'err': 'fqlist_not_found'}
            continue
        fls = get_fq_path(fn_fqlist)
        fls = sorted(fls)
        ct[gse] = len(fls)
        
        if len(fls) > 0:
            for idx, fn in enumerate(fls):
                
                
                
                tmp = f'{gse}\t{fn}\t{kit}\t{seq}'
                res.append(tmp)
                lb = os.path.basename(fn).split('.', 1)[0]
                meta[lb] = [gse, kit, seq]
                if idx < n_fq_per_study:
                    fq_list_for_testing.append(lb)
                    fn_pure = os.path.basename(fn)
                    if node == 'bioinfo2':
                        
                        fn = f'/data/findadapt/download/{fn_pure}'
                    else:
                        fn = f'/data/cqs/chenh19/findadapt/findadapt_bench/download/{fn_pure}'
                    if not os.path.exists(fn):
                        not_found.append(f'{gse}\t{fn}')
                    tmp = f'{gse}\t{fn}\t{kit}\t{seq}'
                    res_lite.append(tmp)

    prefix = 'studies_without_random_seq'
    
    fn_not_found = f'{prefix}.fq_not_found.txt'
    if len(not_found) > 0:
        print(f'{len(not_found)} fastq files not found')
        with open(fn_not_found, 'w') as o:
            print('\n'.join(not_found), file=o)
    else:
        try:
            os.unlink(fn_not_found)
        except:
            pass
    
    
    with open(f'{prefix}.fq_found.txt', 'w') as o:
        print('\n'.join(res), file=o)
    with open(f'{prefix}.fq_found.first{n_fq_per_study}.txt', 'w') as o:
        print('\n'.join(res_lite), file=o)

    with open(f'{prefix}.real_adapter_info.json', 'w') as o:
        json.dump(meta, o, indent=3)

    with open(f'{prefix}.fq_to_run.lb.json', 'w') as o:
        json.dump(fq_list_for_testing, o)

    tmp = json.dumps(ct, indent=3)
    with open(f'{prefix}.count.json', 'w') as o:
        o.write(tmp)


def build_script(fq_list, tool, cmd_template):
    pw_out = f'{tool}_res'
    os.makedirs(pw_out, exist_ok=True)
    os.makedirs(f'{pw_out}/log', exist_ok=True)
    os.makedirs(f'{pw_out}/err', exist_ok=True)
    fn_script_list = f'{tool}.script_list.txt'
    
    print(f'building script for {tool}')
    
    n_total_fq = os.popen(f'cat {fq_list}|wc -l').read().strip()
    with open(fn_script_list, 'w') as o, open(fq_list, 'r') as f:
        idx = 0
        cmds = []
        for i in f:
            fn_raw = i.split('\t', 2)[1].strip()
            fn_fq_pure = fn_raw.rsplit('/', 1)[-1]
            
            if not fn_raw.startswith('/data') or not os.path.exists(fn_raw):
                fn_fq = f'{pw}/download/{fn_fq_pure}'
                if not os.path.exists(fn_fq):
                    print(f'fq not exist: {fn_fq}')
                    continue
            else:
                fn_fq = fn_raw
            fq_lb = fn_fq_pure.replace('.fastq.gz', '').replace('.fastq', '')

            cmd_core = cmd_template.replace('@@fq@@', fn_fq).replace('@@pw@@', pw_out).replace('@@lb@@', fq_lb)
            idx += 1
            cmds.append(f'echo {idx} / {n_total_fq}  - {fq_lb}')
            cmds.append(f'/usr/bin/time -v {cmd_core} > {pw_out}/log/{fq_lb}.log 2> {pw_out}/err/{fq_lb}.err')
            cmds.append('sleep 5')
        print('\n'.join(cmds), file=o)
    print('done')
    return cmds

if __name__ == "__main__":
    
    
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('fq_list', help="""the fq list need to build the script, if not specified, will use the fastq files in studies wihout random sequence""", nargs='?')
    ps.add_argument('lb', help="""the prefix for the fq list, if not specified, will infer from the input name""", nargs='?')
    ps.add_argument('-tool', help="""the software need to build the script, lowercase, nargs=*, default is findadapt only""", choices=['all', 'findadapt', 'dnapi', 'earrings', 'earrings_sensitive', 'atropos_heuristic', 'atropos_kmer'])
    args = ps.parse_args()
    tools_to_run = args.tool or ['findadapt']
    
    if not args.fq_list:
        check_fastq_fls()
        # {'GSE174224': {'not_found': 9, 'found': 0}, 'GSE117744': {'not_found': 2, 'found': 0}, 'GSE123972': {'not_found': 43, 'found': 4}, 'GSE53451': {'not_found': 14, 'found': 0}}
        
        # the access of /nobackup is very slow, so we copy the 362 fastq files to /data
        # cut -f2 2023_realdata/studies_without_random_seq.fq_found.first5.txt|parallel -j5 'cp {} /data/cqs/chenh19/findadapt/findadapt_bench/download'
        
        # real dataset
        fq_list = 'studies_without_random_seq.fq_found.first5.txt'
        prefix = 'studies_without_random_seq'
    else:
        fq_list = args.fq_list
        prefix = args.lb or os.path.basename(fq_list).replace('.fqlist', '').replace('.txt', '')
    
    with open(fq_list) as f:
        lb_list = []
        for i in f:
            fn_raw = i.split('\t', 2)[1].strip()
            fn_fq_pure = fn_raw.rsplit('/', 1)[-1]
            
            if not fn_raw.startswith('/data') or not os.path.exists(fn_raw):
                fn_fq = f'{pw}/download/{fn_fq_pure}'
                fq_lb = fn_fq_pure.split('.', 1)[0]
                if not os.path.exists(fn_fq):
                    print(f'fq not exist: {fn_fq}')
                    continue
            else:
                fn_fq = fn_raw
                fq_lb = fn_fq_pure.replace('.fastq.gz', '').replace('.fastq', '')
            lb_list.append(fq_lb)
    

    with open(f'{prefix}.fq_to_run.lb.json', 'w') as o:
        json.dump(lb_list, o)
    
    
    # @@pw@@ = output pw
    # @@lb@@ = fastq lb
    # @@fq@@ = fastq file full path
    cmd_template = [
        ['findadapt', 'findadapt @@fq@@ -o @@pw@@/@@lb@@'],
        ['DNApi', '/nobackup/h_vangard_1/chenh19/findadapt_bench/tools/DNApi/dnapi.py @@fq@@'],
        ['fastp', 'fastp -i @@fq@@'],
        ['EARRINGS','singularity run /data/cqs/chenh19/dock/earings.sif single -p /nobackup/h_vangard_1/chenh19/findadapt_bench/tools/earrings/hg38 -1 @@fq@@ -o /dev/null'],
        ['EARRINGS_sensitive','singularity run /data/cqs/chenh19/dock/earings.sif single -p /nobackup/h_vangard_1/chenh19/findadapt_bench/tools/earrings/hg38 -1 @@fq@@  --sensitive -o /dev/null'],
        ['atropos_heuristic', 'atropos detect -se @@fq@@ -d heuristic'],
        ['atropos_kmer', 'atropos detect -se @@fq@@ -d khmer'],
    ]
    cmds = []
    for tool, template in cmd_template:
        if 'all' in tools_to_run or tool in tools_to_run:
            cmds += build_script(fq_list, tool, template)

    with open('all.script_list.txt', 'w') as o:
        print('\n'.join(cmds), file=o)
