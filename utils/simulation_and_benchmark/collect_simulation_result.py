"""
collect the result from all tools
1. get the identified adapter, if multiple, choose the one with highest k-mer frequency
2. compare with the read adapter setting, there are 3 category, 3psequence, 5p-random-mer 3p-random-mer, and a final "correct" column, meaning all 3 passed
3. there are 9 rows, which means 9 simulation mode
"""
import os, sys, re
import json


pwout = '/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench'
def get_ground_truth():
    fn = f'{pwout}/simudata/reads_struc.all_setup.txt'
    # 180 lines, 9 settings * 20 datasets
    
    res = {}
    err = {k: [] for k in ['polya not found', 'wrong5p', 'wrong3p', 'wrong empty', 'wrongboth', 'adapter_only']}
    
    polya = 'A' * 12
    with open(fn) as f:
        f.readline()
        for i in f:
            # don't use strip, because the 3padapter at the end can be empty
            tmp = i[:-1].split('\t')
            ifn, len5, len3, adapt = tmp
            len5 = int(len5)
            len3 = int(len3)
            if not adapt:
                adapt = 'empty'
            lb_setup = ifn.split('.')[1]
            covered = 0
            
            if 'polya' in lb_setup:
                covered = 1
                if adapt != "A"*12:
                    err['polya not found'].append(tmp)
            if lb_setup == 'no_adapt':
                covered = 1
                if not (adapt == 'empty' and len5 == 0 and len3 == 0):
                    err['wrong empty'].append(tmp)
            if "add_both_side" in lb_setup:
                covered = 1
                if not (adapt != 'empty' and len5 > 0 and len3 > 0):
                    err['wrongboth'].append(tmp)
                    

            if 'adapter_only' in lb_setup:
                covered = 1
                if not (adapt != 'empty' and adapt != polya  and len5 == 0 and len3 == 0):
                    err['adapter_only'].append(tmp)
                    
            if 'add_5mer' in lb_setup:
                covered = 1
                if not (adapt != 'empty' and len5 > 0 and len3 == 0):
                    err['wrong5p'].append(tmp)

            if 'add_3mer' in lb_setup:
                covered = 1
                if not (adapt != 'empty' and len5 == 0 and len3 > 0):
                    err['wrong3p'].append(tmp)


            if not covered:
                print(f'not covered: {tmp}')
            
            
            lb_data = os.path.basename(ifn).replace('.fastq', '')
            res[lb_data] = [lb_setup, len5, len3, adapt]
    
    # print(json.dumps(res, indent=3))
    
    err = {k: v for k, v in err.items() if len(v) > 0}
    if len(err) > 0:
        print(f'the following dataset has problem')
        print([(k, len(v)) for k, v in err.items()])
        with open('{pwout}/bench_res/datasets_with_problem.err', 'w') as o:
            for k, v in err.items():
                for line in v:
                    print('\t'.join([k] + line), file=o)
    return res

def process_atropos(data_lb, fn, setup):
    lb, len5_exp, len3_exp, adapt_exp = setup
    # k_tool = ='atropos_kmer'  or atropos_heuristic
    
    
    with open(fn) as f:
        content = f.read()
        # Detected 0 adapters/contaminants:
        # Detected 20 adapters/contaminants:
        # 1. Longest kmer: AACTGTAGGCAC
        #     Abundance (full-length) in 10000 reads: 5864 (58.6%)
        #     Frequency of k-mers: 2.35%
        len5_real = 0
        len3_real = 0
        adapt_real = None
        
        f.seek(0)
        if "Detected 0 adapters/contaminants" in content:
            adapt_real = 'empty'
        else:
            for line in f:
                if 'Longest kmer' in line:
                    adapt_real = line.rsplit(':', 1)[-1].strip()
                    break

        
    if adapt_real is None:
        print(fn)
        return 0, 'NA'
    if adapt_real == adapt_exp and len5_real == len5_exp and len3_real == len3_exp:
        matched = 1
    else:
        matched = 0
    matched_str = 'match' if matched else 'unmatch'
    detail = [matched_str, f'{len5_exp}:{len5_real}', f'{len3_exp}:{len3_real}', f'{adapt_exp}:{adapt_real}', data_lb]
    
    return matched, detail


def process_findadapt(data_lb, fn, setup):
    lb, len5_exp, len3_exp, adapt_exp = setup
    # k_tool = ='atropos_kmer'  or atropos_heuristic
    with open(fn) as f:
        f.readline()
        line = f.readline().strip().split('\t')
        adapt_real = line[2]
        len3_real = int(line[3])
        len5_real = int(line[6])
        if not adapt_real:
            adapt_real = 'empty'

    if adapt_real == adapt_exp and len5_real == len5_exp and len3_real == len3_exp:
        matched = 1
    else:
        matched = 0
    matched_str = 'match' if matched else 'unmatch'
    detail = [matched_str, f'{len5_exp}:{len5_real}', f'{len3_exp}:{len3_real}', f'{adapt_exp}:{adapt_real}', data_lb]
    
    return matched, detail



def process_earrings(data_lb, fn, setup):
    lb, len5_exp, len3_exp, adapt_exp = setup
    # k_tool = ='atropos_kmer'  or atropos_heuristic
    with open(fn) as f:
        content = f.read()
        len5_real = 0
        len3_real = 0
        f.seek(0)
        if "unable to detect adapter" in content:
            adapt_real = 'empty'
        else:
            collect = 0
            for line in f:
                # adapter found: GAAAAAAAAAAAAA
                if 'adapter found' in line:
                    adapt_real = line.rsplit(':', 1)[-1].strip()

    if adapt_real == adapt_exp and len5_real == len5_exp and len3_real == len3_exp:
        matched = 1
    else:
        matched = 0
    matched_str = 'match' if matched else 'unmatch'
    detail = [matched_str, f'{len5_exp}:{len5_real}', f'{len3_exp}:{len3_real}', f'{adapt_exp}:{adapt_real}', data_lb]
    
    return matched, detail


def process_fastp(data_lb, fn, setup):
    lb, len5_exp, len3_exp, adapt_exp = setup
    # k_tool = ='atropos_kmer'  or atropos_heuristic
    with open(fn) as f:
        content = f.read()
        # Detected 0 adapters/contaminants:
        # Detected 20 adapters/contaminants:
        # 1. Longest kmer: AACTGTAGGCAC
        #     Abundance (full-length) in 10000 reads: 5864 (58.6%)
        #     Frequency of k-mers: 2.35%
        len5_real = 0
        len3_real = 0
        f.seek(0)
        if "No adapter detected" in content:
            adapt_real = 'empty'
        else:
            collect = 0
            for line in f:
                if collect == 0 and 'Detecting adapter' in line:
                    collect = 1
                    continue
                if collect:
                    tmp = line.strip()
                    if not re.match('[ATCG]+$', tmp):
                        print(f'ERROR, inlalid format for {fn}')
                        adapt_real = 'NA'
                    else:
                        adapt_real = tmp
                    break

    if adapt_real == adapt_exp and len5_real == len5_exp and len3_real == len3_exp:
        matched = 1
    else:
        matched = 0
    matched_str = 'match' if matched else 'unmatch'
    detail = [matched_str, f'{len5_exp}:{len5_real}', f'{len3_exp}:{len3_real}', f'{adapt_exp}:{adapt_real}', data_lb]
    
    return matched, detail


def process_dnapi(data_lb, fn, setup):
    lb, len5_exp, len3_exp, adapt_exp = setup
    # k_tool = ='atropos_kmer'  or atropos_heuristic
    with open(fn) as f:
        len5_real = 0
        len3_real = 0
        adapt_real = f.readline().strip()
        

    if adapt_real == adapt_exp and len5_real == len5_exp and len3_real == len3_exp:
        matched = 1
    else:
        matched = 0
    matched_str = 'match' if matched else 'unmatch'
    detail = [matched_str, f'{len5_exp}:{len5_real}', f'{len3_exp}:{len3_real}', f'{adapt_exp}:{adapt_real}', data_lb]
    
    return matched, detail



def main():
    true_res = get_ground_truth()
    
    n = 0

    setup_lb = set([_[0] for _ in true_res.values()])
    tools = [
            'findadapt',
            'atropos_kmer',
            'atropos_heuristic',
            'earrings_hg38',
            'earrings_hg38_sensitive',
            'fastp',
            'dnapi',
    ]
    
    
    process_func_map = {
            'findadapt': process_findadapt,
            'atropos_kmer': process_atropos,
            'atropos_heuristic': process_atropos,
            'earrings_hg38': process_earrings,
            'earrings_hg38_sensitive': process_earrings,
            'fastp': process_fastp,
            'dnapi': process_dnapi,
    }
    
    # tools = ['earrings_hg38', 'earrings_hg38_sensitive']
    
    detail = {k: [] for k in tools}
    res = {k: {k1: [] for k1 in tools} for k in setup_lb}  # k1 = setup label, k2 = tool, v = list of 0 / 1, if match with setup, will be 1, else 0
    # print(setup_lb)
    # {'polya_add_5mer', 'polya_only', 'no_adapt', 'polya_add_both_side', 'adapter_only', 'polya_add_3mer', 'add_5mer', 'add_3mer', 'add_both_side'}
    
    res_not_found = []
    for data_lb, setup in true_res.items():
        # if n > 1:
        #     break
        # setup = [lb_setup, len5, len3, adapt]

        n += 1
        lb = setup[0]
        # setup[1] = str(setup[1])
        # setup[2] = str(setup[2])
        
        for k_tool in tools:
            ires = res[lb][k_tool]
            idetail = detail[k_tool]
            process_func = process_func_map[k_tool]
            
            if k_tool == 'findadapt':
                fn = f'{pwout}/bench_res/{k_tool}/{data_lb}.adapter.txt'
            else:
                fn = f'{pwout}/bench_res/{k_tool}/{data_lb}.log'
            if not os.path.exists(fn):
                res_not_found.append([data_lb, k_tool, fn])
                continue 
            
            matched, detail_line = process_func(data_lb, fn, setup)
            ires.append(matched)
            idetail.append(detail_line)
    
    # print(res)
    
    if len(res_not_found) > 0:
        print(f'the following result files are not found: \n' + '\n'.join(map(str, res_not_found)))

    for toolname, lines in detail.items():
        lines = sorted(lines, key=lambda _: (_[0], _[4]))
        lines = ['\t'.join(_) for _ in lines]
        fn_detail = f'{pwout}/bench_res/detail.{toolname}.txt'
        with open(fn_detail,'w') as o:
            print('\n'.join(lines), file=o)

    with open(f'{pwout}/bench_res/benchmark_final_res.json', 'w') as o:
        json.dump(res, o, indent=3)

main()
