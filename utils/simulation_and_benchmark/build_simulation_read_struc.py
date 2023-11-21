#! /usr/bin/env python3
"""
build the simulation data with different read structure
1) simple, only 3' adapter 
2) 1-5bp random-mer at 5' end
3) 1-5bp random-mer at 3' end
4) 1-5bp random-mer at both end
5) ploly(A) as 3' adapter
6) 1-5bp random-mer at 5'end and poly(A) tail
7) 1-5bp random-mer at 3' end and poly(A)
8) 1-5bp randomer at both end and poly(A)
9) no adapters

20 datasets for each scenario

"""

import os, sys, re
import gzip
import json
import pickle
from random import choice, sample, random
bases = 'ATCG'
random_len_range = range(1, 6)
adapt_pool = [
    'AGATCGGAAGAG',  # NEBNext
    'AACTGTAGGCAC', # kiaseq
    'TGGAATTCTCGG', # nextflex, truseq
]

other_bases = {base: bases.replace(base, '') for base in bases}

def add_mutation(adapt_seq, mut_ratio=0.0025):
    # the mutation ratio = 0.25%
    # each base hs the possibility of mutate
    adapt_seq_new = []
    n_mut = 0
    for base in adapt_seq:

        point = random()
        if point < mut_ratio:
            n_mut += 1
            adapt_seq_new.append(choice(other_bases[base]))
        else:
            adapt_seq_new.append(base)

    return ''.join(adapt_seq_new), n_mut

def load_posttrim_data():
    # GSE109356
    # /nobackup/h_vangard_1/chenh19/exRNA/test/download/SRR6502962/SRR6502962.fastq.gz
    
    # cutadapt -j4 -q 20 -a AGATCGGAAGAG -O 1 --trim-n -m 16 -o SRR6502962.trimmed.fastq.gz /gpfs23/scratch/h_vangard_1/chenh19/exRNA/test/download/SRR6502962/SRR6502962.fastq.gz
    
        # Total reads processed:                 963,413
        # Reads with adapters:                   958,326 (99.5%)

        # == Read fate breakdown ==
        # Reads that were too short:              36,330 (3.8%)
        # Reads written (passing filters):       927,083 (96.2%)

        # Total basepairs processed:    48,170,650 bp
        # Quality-trimmed:                  57,405 bp (0.1%)
        # Total written (filtered):     25,703,471 bp (53.4%)

    # 15M , 1,023,493 reads
    # fn = '/data/cqs/chenh19/findadapt/findadapt_bench/simu_read_structure/SRR6502962.trimmed.fastq.gz'
    fn = '/data/cqs/chenh19/findadapt/findadapt_bench/simu_read_structure/SRR6502962.trimmed.fastq.gz'
    
    fn_pkl = fn.replace('.fastq.gz', '.fastq.pkl')
    
    if os.path.exists(fn_pkl):
        with open(fn_pkl, 'rb') as f:
            return pickle.load(f)
    
    print('parsing the fastq to python dict')
    res = {}
    with gzip.open(fn, 'rt') as f:
        n = 0
        while True:
            name = f.readline().strip()
            if not name:
                break
            seq = f.readline().strip()
            f.readline()
            qscore = f.readline().strip()
            res[n] = [seq, qscore, name]
            n += 1
    print(f'parsed reads = {n}')
    with open(fn_pkl, 'wb') as o:
        pickle.dump(res, o)
    return res


rawreads = load_posttrim_data()
sn_all = range(len(rawreads))


def get_random_mer(length):
    return ''.join([choice(bases) for _ in range(length)])


def get_write_read_func(setup):
    len_5p, len_3p, adapt_seq = setup
    len_adapt = len(adapt_seq)
    bases = 'ATCG'
    reads_total_len = 50
    mut_ratio = 0.0025
    
    def fill_read_with_random(seq):
        """
        fill the reads with random base to make the final reads length as reads_total_len
        """
        n_fill = reads_total_len - len_5p - len_3p - len_adapt - len(seq)
        if n_fill == 0:
            return '', ''
        return ''.join([choice(bases) for _ in range(n_fill)]), 'F' * n_fill
        
    
    
    if len_5p == 0 and len_3p == 0 and len_adapt > 0:
        def inner_func(iread):
            seq, qscore, header = iread
            adapt_seq_new, n_mut = add_mutation(adapt_seq, mut_ratio)
            seq_fill, qscore_fill = fill_read_with_random(seq)
            seqnew = seq + adapt_seq_new + seq_fill
            qscore_new = qscore + 'F' * len_adapt + qscore_fill
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires), n_mut
        return inner_func

    if len_5p > 0 and len_3p == 0 and len_adapt > 0:
        def inner_func(iread):
            seq, qscore, header = iread
            seq_5p = ''.join([choice(bases) for _ in range(len_5p)])
            adapt_seq_new, n_mut = add_mutation(adapt_seq, mut_ratio)
            seq_fill, qscore_fill = fill_read_with_random(seq)
            seqnew = seq_5p + seq + adapt_seq_new + seq_fill
            qscore_new = 'F' * len_5p + qscore + 'F' * len_adapt + qscore_fill
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires), n_mut
        return inner_func
    

    if len_5p == 0 and len_3p > 0 and len_adapt > 0:
        len1 = len_3p + len_adapt
        def inner_func(iread):
            seq, qscore, header = iread
            seq_fill, qscore_fill = fill_read_with_random(seq)
            seq_3p = ''.join([choice(bases) for _ in range(len_3p)])
            adapt_seq_new, n_mut = add_mutation(adapt_seq, mut_ratio)
            seqnew = seq + seq_3p + adapt_seq_new + seq_fill
            qscore_new = qscore + 'F' * len1 + qscore_fill
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires), n_mut
        return inner_func
    

    if len_5p > 0 and len_3p > 0 and len_adapt > 0:
        len1 = len_3p + len_adapt
        def inner_func(iread):
            seq, qscore, header = iread
            seq_fill, qscore_fill = fill_read_with_random(seq)
            seq_5p = ''.join([choice(bases) for _ in range(len_5p)])
            seq_3p = ''.join([choice(bases) for _ in range(len_3p)])
            adapt_seq_new, n_mut = add_mutation(adapt_seq, mut_ratio)
            seqnew = seq_5p + seq + seq_3p + adapt_seq_new + seq_fill
            qscore_new = 'F' * len_5p + qscore + 'F' * len1 + qscore_fill
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires), n_mut
        return inner_func
    

    if len_5p == 0 and len_3p == 0 and len_adapt == 0:
        def inner_func(iread):
            seq, qscore, header = iread
            ires = [header, seq, '+', qscore]
            return '\n'.join(ires), 0
        return inner_func


def write_file(fno, setup):
    # rawreads is the dict got from load_posttrim_data
    # key = number, the sn of the reads
    # v = list  [post_trim_seq,  qscore, read_header]
    # the read header has the newline symbol
    
    func_write_read = get_write_read_func(setup)
    
    len_5p, len_3p, adapt_seq = setup
    n_reads_out = 10_000

    sn_picked = sample(sn_all, n_reads_out)
    
    fnconfig = fno.replace('.fastq', '.setup.txt')
    with open(fnconfig, 'w') as o:
        print('5p_phase\t3p_phase\t3p_adapter', file=o)
        print(f'{len_5p}\t{len_3p}\t{adapt_seq}', file=o)
    
    last_base_freq = {_: 0 for _ in 'ATCG'}
    n_mut_in_adapter_all = {}
    
    
    with open(fno, 'w') as o:
        for sn in sn_picked:
            # ['TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT',
            # '?????B<B?9BDDD?BFFFFFFBFHHHH?CCHH',
            # '@SRR6502962.1 1/1\n']
            last_base = rawreads[sn][0][0]
            last_base_freq[last_base] += 1
            new_read, n_mut_in_adapter = func_write_read(rawreads[sn])
            print(new_read, file=o)
            n_mut_in_adapter_all.setdefault(n_mut_in_adapter, 0)
            n_mut_in_adapter_all[n_mut_in_adapter] += 1
    
    
    return last_base_freq, n_mut_in_adapter_all

last_base_freq_all = {}

def main(iter_sn, add_5mer, add_3mer, add_polya, add_3p_adapt):
    # 
    if add_5mer:
        length_5p = choice(random_len_range)
    else:
        length_5p = 0
    
    if add_3mer:
        length_3p = choice(random_len_range)
    else:
        length_3p = 0
    

    if add_polya:
        adapt_seq = 'A' * 12
    elif add_3p_adapt:
        adapt_seq = choice(adapt_pool)
    else:
        adapt_seq = ''
    
    # record the setup
    
    # 1) simple, only 3' adapter 
    # 2) 1-5bp random-mer at 5' end
    # 3) 1-5bp random-mer at 3' end
    # 4) 1-5bp random-mer at both end
    # 5) ploly(A) as 3' adapter
    # 6) 1-5bp random-mer at 5'end and poly(A) tail
    # 7) 1-5bp random-mer at 3' end and poly(A)
    # 8) 1-5bp randomer at both end and poly(A)
    # 9) no adapters

    if add_5mer == False and add_3mer == False and add_3p_adapt == True:
        read_struc_type = 'adapter_only'
    elif add_5mer == True and add_3mer == False and add_3p_adapt == True:
        read_struc_type = 'add_5mer'
    elif add_5mer == False and add_3mer == True and add_3p_adapt == True:
        read_struc_type = 'add_3mer'
    elif add_5mer == True and add_3mer == True and add_3p_adapt == True:
        read_struc_type = 'add_both_side'
    elif add_5mer == False and add_3mer == False and add_polya == True:
        read_struc_type = 'polya_only'
    elif add_5mer == True and add_3mer == False and add_polya == True:
        read_struc_type = 'polya_add_5mer'
    elif add_5mer == False and add_3mer == True and add_polya == True:
        read_struc_type = 'polya_add_3mer'
    elif add_5mer == True and add_3mer == True and add_polya == True:
        read_struc_type = 'polya_add_both_side'
    elif add_5mer == False and add_3mer == False and add_polya == False and add_3p_adapt == False:
        read_struc_type = 'no_adapt'
    
    
    # fix the previous error
    # modify here
    # if 'add_5mer' not in read_struc_type and 'add_3mer' not in read_struc_type:
    #     return

    
    lb = f'{read_struc_type}.{iter_sn}'
    if iter_sn == 1:
        print(read_struc_type)
    fno = f'{pwout}/simu.{lb}.fastq'
    simu_setup = [length_5p, length_3p, adapt_seq]
    
    last_base_freq, n_mut_in_adapter = write_file(fno, simu_setup)
    
    last_base_freq_all[lb] = last_base_freq
    n_mut_in_adapter_by_file[lb] = n_mut_in_adapter

    tmp = [read_struc_type, read_struc_type, adapt_seq]

    return '\t'.join([read_struc_type, fno, str(length_5p), str(length_3p), adapt_seq]), f'simu.{lb}', tmp

if __name__ == "__main__":
    
    # pwout = '/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/simudata/reads_struc'
    pwout = '/data/cqs/chenh19/findadapt/findadapt_bench/simu_read_structure/reads_struc'
    
    iter_times = 20
    
    setup_list = []
    n_mut_in_adapter_by_file = {}
    
    real_adapter_info = {}
    
    for iter_sn in range(iter_times):
        # noadapter
        ires, simu_lb, tmp = main(iter_sn+1, False, False, False, False)
        real_adapter_info[simu_lb] = tmp
        if ires:
            setup_list.append(ires)

    for adapt_seq in [[True, False], [False, True]]:
        # add_polya, add_3p_adap
        for add_5mer in [True, False]:
            for add_3mer in [True, False]:
                for iter_sn in range(iter_times):
                    ires, simu_lb, tmp = main(iter_sn+1, add_5mer, add_3mer, *adapt_seq)
                    real_adapter_info[simu_lb] = tmp
                    if ires:
                        setup_list.append(ires)

    # save the all simulation setup
    
    # with open('/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/simudata/reads_struc.all_setup.txt', 'w') as o:
        # print('fn\tlen5\tlen3\tadapter', file=o)
    
    with open('/data/cqs/chenh19/findadapt/findadapt_bench/simu_read_structure/simu_struc.real_adapter_info.json', 'w') as o:
        json.dump(real_adapter_info,  o, indent=3)
    
    
    # modify here
    # with open('/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/simudata/reads_struc.all_setup.txt', 'a') as o:
    with open('/data/cqs/chenh19/findadapt/findadapt_bench/simu_read_structure/reads_struc.all_setup.txt', 'w') as o:
        print('\n'.join(setup_list), file=o)
    

    tmp = json.dumps(n_mut_in_adapter_by_file, indent=4)
    print(tmp)
    with open('/data/cqs/chenh19/findadapt/findadapt_bench/simu_read_structure/simulation_read_struc.mut_adapter_count.json', 'w') as o:
        print(tmp, file=o)

    with open('/data/cqs/chenh19/findadapt/findadapt_bench/simu_read_structure/simu_data.insert_seq.last_base_freq.txt', 'w') as o:
        for k in sorted(last_base_freq_all):
            v = last_base_freq_all[k]
            v_sum = sum(v.values())
            v_ratio = {base: round(ct/v_sum, 3) for base, ct in v.items()}
            print(f'{k}\t{v}\t{v_ratio}', file=o)

