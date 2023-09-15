#! /usr/bin/env python3
"""
used to collect the run time and memory usage from /usr/bin/time -v
"""

import os, sys, re
import json
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('pw', help="""path for the results, should contain the log and err subfolder""")
ps.add_argument('-rm', help="""remove the result not in expected file list""", action='store_true')
ps.add_argument('-prefix', '-p',  help="""the prefix of the current analysis, if not specified, will infer from the json file name""")
args = ps.parse_args()
import platform
node = platform.node()
if 'bioinfo2' in node:
    pw = '/data/findadapt'
elif 'hanuman' in node:
    pw = '/data/cqs/chenh19/findadapt/findadapt_bench'

red = '\u001b[31;1m'
green = '\u001b[32;1m'
reset = '\u001b[0m'

def parse_time(fn, fq_for_testing):
    with open(fn) as f:
        lb = os.path.basename(fn).replace('.fastq.gz', '').replace('.fastq', '').replace('.err', '').replace('.log', '')
        if lb not in fq_for_testing:
            print(f'file not in test list: {fn}, lb = {lb}')
            
            return None
        cmd = f.readline()
        fq_path = f'{pw}/download/{lb}.fastq.gz'
        
        ires = {'time': None, 'memory': None}
        for i in f:
            i = i.strip()
            if 'wall clock' in i:
                runtime_raw = i.split(': ', 1)[-1].strip()
                runtime_raw = [float(_) for _ in runtime_raw.split(':')]
                runtime = 0
                for stage, num in enumerate(runtime_raw[::-1]):
                    runtime += num * 60**stage

                ires['time'] = runtime
            elif i.startswith('Maximum resident set'):
                mem = float(i.split(':')[-1].strip()) / 1024
                ires['memory'] = mem
    
    return [lb, fq_path, ires['time'], ires['memory']] + real_adapter_info[lb]


def parse_findadapt(pwo, lb):
    fn = f'{pwo}/{lb}.adapter.txt'
    seq_exp = real_adapter_info[lb][2].split(' ')  # gse, kit_exp, seq_3p_exp
    with open(fn) as f:
        # prj	total_reads	3p_seq	3p_phase	3p_count	3p_ratio	5p_phase	5p_count	5p_ratio	err
        f.readline()
        line = f.readline()[:-1].split('\t')
        # seq_3p, phase_3p, ct_3p, phase_5p, ct_5p
        ires = [line[_] for _ in [2, 3, 4, 6, 7]]
        if not ires[0]:
            ires[0] = 'empty'
        if ires[0] in seq_exp and ires[1] == '0' and ires[3] == '0':
            ires.append(f'{green}match{reset}')
        else:
            ires.append(f'{red}fail{reset}')
        
        return ires
    

def parse_dnapi(pwo, lb):
    # DNApi_res/log/SRX8839335.log
    # the content is a single seq
    # the length is not fixed
    fn = f'{pwo}/log/{lb}.log'
    seq_exp = real_adapter_info[lb][2].split(' ')  # gse, kit_exp, seq_3p_exp
    
    with open(fn) as f:
        seq = f.read().strip().replace('\n', '@')
        seq_first_12 = seq[:12]
        ires = [seq, 0, 'NA', 0, 'NA']
        if seq_first_12 in seq_exp or seq[:10] == 'A' * 10:
            ires.append(f'{green}match{reset}')
        else:
            ires.append(f'{red}fail{reset}')
        return ires

def parse_fastp(pwo, lb):

    # Detecting adapter sequence for read1...
    # >Illumina TruSeq Adapter Read 1
    # AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    
    # Detecting adapter sequence for read1...
    # No adapter detected for read1
    
    fn = f'{pwo}/err/{lb}.err'
    seq_exp = real_adapter_info[lb][2].split(' ')  # gse, kit_exp, seq_3p_exp

    with open(fn) as f:
        collect = 0
        for i in f:
            i = i.strip()
            if i.startswith('Detecting adapter sequence'):
                collect = 1
                continue
            elif collect == 1:
                if not i:
                    break
                if i.startswith('No adapter detected'):
                    seq = 'empty'
                    break
                elif i[0] == '>':
                    continue
                else:
                    seq = i
                    break
        
        ires = [seq, 0, 'NA', 0, 'NA']
        seq_first_12 = seq[:12]
        if seq_first_12 in seq_exp or seq[:10] == 'A' * 10:
            ires.append(f'{green}match{reset}')
        else:
            ires.append(f'{red}fail{reset}')
        return ires
    
def parse_earring(pwo, lb):
    # Index prefix: /nobackup/h_vangard_1/chenh19/findadapt_bench/tools/earrings/hg38
    # Input file name: /data/cqs/chenh19/findadapt/findadapt_bench/download/SRR10196215.fastq.gz
    # Output file name: /dev/null.fastq
    # # of threads: 1
    # Is fastq: true, Is gz input: true, Is bam: false
    # Seed length: 50, Max alignment: 0, No mismatch: false
    # Prune factor: 0.03, Sensitive mode: false
    # Min length: 0, UMI: false
    # Default adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    # load table success
    # unable to detect adapter, use default adapter
    # 10.3442 sec
    
    # load table success
    # adapter found: TCATTTCTTTGTT
    # 8.7337 sec
    fn = f'{pwo}/log/{lb}.log'
    seq_exp = real_adapter_info[lb][2].split(' ')  # gse, kit_exp, seq_3p_exp

    with open(fn) as f:
        for i in f:
            i = i.strip()
            if i.startswith('Default adapter:'):
                default_seq = i.split(':')[1].strip()
            elif i.startswith('unable to detect adapter'):
                seq_raw = 'default'
            elif i.startswith('adapter found'):
                seq_raw = i.split(':')[1].strip()
        
        if seq_raw == 'default':
            seq = default_seq + " (default)"
        else:
            seq = seq_raw
        
        ires = [seq, 0, 'NA', 0, 'NA']
        seq_first_12 = seq[:12]
        if seq_first_12 in seq_exp or seq[:10] == 'A' * 10:
            ires.append(f'{green}match{reset}')
        else:
            ires.append(f'{red}fail{reset}')
        return ires

    #  the fastq files that EARRINGS found the adapter
    # ls EARRINGS_res/log/*|while read fn;do t=$(grep "unable to detect adapter" $fn); if [[ -z "$t" ]];then echo $fn;fi;done|wc -l
    

def parse_atropos(pwo, lb):
    # =======
    # Input 1
    # =======

    # File: /data/cqs/chenh19/findadapt/findadapt_bench/download/SRR10196216.fastq.gz
    # Detected 2 adapters/contaminants:
    # 1. Longest kmer: CTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGCTAGCTGC
    # Longest matching sequence: CTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGCTAGCTGC
    # Abundance (full-length) in 10000 reads: 33 (0.3%)
    # Number of k-mer matches: 53086
    # 2. Longest kmer: GTATAGTGGTGAGTATCCCCGCCTGTCACGCGGGAGACCGGGGTTCGATTCC
    # Longest matching sequence: GTATAGTGGTGAGTATCCCCGCCTGTCACGCGGGAGACCGGGGTTCGATTCCC
    # Abundance (full-length) in 10000 reads: 24 (0.2%)
    # Number of k-mer matches: 33481


    # k-mer mode
    # =======
    # Input 1
    # =======

    # File: /data/cqs/chenh19/findadapt/findadapt_bench/download/SRR10196216.fastq.gz
    # Detected 3 adapters/contaminants:
    # 1. Longest kmer: CTACCACTGAAC
    # Frequency of k-mers: 0.14%
    # 2. Longest kmer: CCACTGAACCAC
    # Frequency of k-mers: 0.12%
    # 3. Longest kmer: ATTCTACCACTG
    # Frequency of k-mers: 0.12%

    fn = f'{pwo}/log/{lb}.log'
    seq_exp = real_adapter_info[lb][2].split(' ')  # gse, kit_exp, seq_3p_exp

    with open(fn) as f:
        seq = 'empty'
        for i in f:
            i = i.strip()
            if 'Longest kmer' in i:
                seq = i.split(':')[-1].strip()
                break
        ires = [seq, 0, 'NA', 0, 'NA']
        seq_first_12 = seq[:12]
        if seq_first_12 in seq_exp or seq[:10] == 'A' * 10:
            ires.append(f'{green}match{reset}')
        else:
            ires.append(f'{red}fail{reset}')
        return ires


if __name__ == "__main__":
    pwo = args.pw
    pw_err = f'{pwo}/err'
    fls = os.listdir(f'{pw_err}')
    rm_extra = args.rm
    
    tool = os.path.basename(pwo).replace('_res', '')
    
    res = []
    
    prefix = args.prefix
    if not prefix:
        from glob import glob
        tmp = glob('*.fq_to_run.lb.json')
        if len(tmp) == 0:
            print('ERROR, .fq_to_run.lb.json not exist')
            sys.exit(1)
        prefix = tmp[0].replace('.fq_to_run.lb.json', '')

    with open(f'{prefix}.fq_to_run.lb.json') as f:
        fq_for_testing = json.load(f)
        fq_for_testing = set(fq_for_testing)
    
    fn_real_adapter_info = f'{prefix}.real_adapter_info.json'
    if os.path.exists(fn_real_adapter_info):
        with open(fn_real_adapter_info) as f:
            real_adapter_info = json.load(f)
    else:
        real_adapter_info = {}
        for lb in fq_for_testing:
            # gse, kit_exp, seq_3p_exp
            real_adapter_info[lb] = ['NA', 'NA', 'NA']


    func_parse_adapter_res = {
        'findadapt':  parse_findadapt,
        'dnapi': parse_dnapi,
        'fastp': parse_fastp,
        'earrings': parse_earring,
        'earrings_sensitive': parse_earring,
        'atropos_kmer': parse_atropos,
        'atropos_heuristic': parse_atropos,
    }[tool.lower()]
    
    for fn in fls:
        fn = f'{pw_err}/{fn}'
        res_time_mem = parse_time(fn, fq_for_testing)
        if not res_time_mem:
            if rm_extra:
                os.unlink(fn)
            continue
        
        lb = res_time_mem[0]
        # get the result
        
        res_adapter = func_parse_adapter_res(pwo, lb)
        if res_adapter is None:
            print(f'ERROR, fail to get adapter for {lb}')
            continue
        
        res.append(res_time_mem + res_adapter)

    header = 'lb,fq,wall_time,max_mem,gse,kit_exp,seq_3p_exp,seq_3p,phase_3p,ct_3p,phase_5p,ct_5p,match'.replace(',', '\t')
    
    res = sorted(res, key=lambda _: (_[4], _[0]))
    prj_lb = os.path.basename(pwo)
    with open(f'{prj_lb}.time_mem_summary.txt', 'w') as o:
        print(header, file=o)
        tmp = ['\t'.join(map(str, _)) for _ in res]
        print('\n'.join(tmp), file=o)

    
    