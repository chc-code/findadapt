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
import pickle
from random import choice, sample
bases = 'ATCG'
random_len_range = range(1, 6)
adapt_pool = [
    'AGATCGGAAGAG',  # NEBNext
    'AACTGTAGGCAC', # kiaseq
    'TGGAATTCTCGG', # nextflex, truseq
]

def load_posttrim_data():
    # GSE109356
    # /gpfs23/scratch/h_vangard_1/chenh19/exRNA/test/download/SRR6502962/SRR6502962.fastq.gz
    
    # cutadapt -j4 -q 20 -a AGATCGGAAGAG -O 1 --trim-n -m 16 -o SRR6502962.trimmed.fastq.gz /gpfs23/scratch/h_vangard_1/chenh19/exRNA/test/download/SRR6502962/SRR6502962.fastq.gz
    
        # Total reads processed:                 963,413
        # Reads with adapters:                   958,326 (99.5%)

        # == Read fate breakdown ==
        # Reads that were too short:              36,330 (3.8%)
        # Reads written (passing filters):       927,083 (96.2%)

        # Total basepairs processed:    48,170,650 bp
        # Quality-trimmed:                  57,405 bp (0.1%)
        # Total written (filtered):     25,703,471 bp (53.4%)

    fn = '/scratch/h_vangard_1/chenh19/findadapt_bench/simudata/SRR6502962.trimmed.fastq.gz'
    
    fn_pkl = fn.replace('.fastq.gz', '.fastq.pkl')
    
    if os.path.exists(fn_pkl):
        with open(fn_pkl, 'rb') as f:
            return pickle.load(f)
    
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
    if len_5p == 0 and len_3p == 0 and len_adapt > 0:
        def inner_func(iread):
            seq, qscore, header = iread
            seqnew = seq + adapt_seq
            qscore_new = qscore + 'F' * len_adapt
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires)
        return inner_func

    if len_5p > 0 and len_3p == 0 and len_adapt > 0:
        def inner_func(iread):
            seq, qscore, header = iread
            seq_5p = ''.join([choice(bases) for _ in range(len_5p)])
            seqnew = seq_5p + seq + adapt_seq
            qscore_new = 'F' * len_5p + qscore + 'F' * len_adapt
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires)
        return inner_func
    

    if len_5p == 0 and len_3p > 0 and len_adapt > 0:
        len1 = len_3p + len_adapt
        def inner_func(iread):
            seq, qscore, header = iread
            seq_3p = ''.join([choice(bases) for _ in range(len_3p)])
            seqnew = seq + seq_3p + adapt_seq
            qscore_new = qscore + 'F' * len1
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires)
        return inner_func
    

    if len_5p > 0 and len_3p > 0 and len_adapt > 0:
        len1 = len_3p + len_adapt
        def inner_func(iread):
            seq, qscore, header = iread
            seq_5p = ''.join([choice(bases) for _ in range(len_5p)])
            seq_3p = ''.join([choice(bases) for _ in range(len_3p)])
            seqnew = seq_5p + seq + seq_3p + adapt_seq
            qscore_new = 'F' * len_5p + qscore + 'F' * len1
            ires = [header, seqnew, '+', qscore_new]
            return '\n'.join(ires)
        return inner_func
    

    if len_5p == 0 and len_3p == 0 and len_adapt == 0:
        def inner_func(iread):
            seq, qscore, header = iread
            ires = [header, seq, '+', qscore]
            return '\n'.join(ires)
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
    
    with open(fno, 'w') as o:
        for sn in sn_picked:
            # ['TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT',
            # '?????B<B?9BDDD?BFFFFFFBFHHHH?CCHH',
            # '@SRR6502962.1 1/1\n']
            print(func_write_read(rawreads[sn]), file=o)
            


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
    if 'add_5mer' not in read_struc_type and 'add_3mer' not in read_struc_type:
        return
    
    
    
    lb = f'{read_struc_type}.{iter_sn}'
    if iter_sn == 1:
        print(read_struc_type)
    fno = f'{pwout}/simu.{lb}.fastq'
    simu_setup = [length_5p, length_3p, adapt_seq]
    write_file(fno, simu_setup)
    
    return '\t'.join([fno, str(length_5p), str(length_3p), adapt_seq])

if __name__ == "__main__":
    
    pwout = '/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/simudata/reads_struc'
    
    iter_times = 20
    
    setup_list = []
    for iter_sn in range(iter_times):
        # noadapter
        ires = main(iter_sn+1, False, False, False, False)
        if ires:
            setup_list.append(ires)

    for adapt_seq in [[True, False], [False, True]]:
        # add_polya, add_3p_adap
        for add_5mer in [True, False]:
            for add_3mer in [True, False]:
                for iter_sn in range(iter_times):
                    ires = main(iter_sn+1, add_5mer, add_3mer, *adapt_seq)
                    if ires:
                        setup_list.append(ires)

    
    # save the all simulation setup
    
    # with open('/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/simudata/reads_struc.all_setup.txt', 'w') as o:
        # print('fn\tlen5\tlen3\tadapter', file=o)
    
    # modify here
    with open('/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/simudata/reads_struc.all_setup.txt', 'a') as o:
        
        
        print('\n'.join(setup_list), file=o)
    
    
    
