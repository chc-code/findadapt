import os
def get_gse_list():
    fn = '/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/realdata/gse_list.87_datasets.txt'
    
    with open(fn) as f:
        return [_.strip() for _ in f if _.strip()]

def get_fq(fn):
    if os.path.exists(fn):
        return fn
    
    for ipw in ['geo_pe', 'geo_se', 'liu_xiao']:
        ifn = fn.replace('/test/', f'/{ipw}/')
        
        if os.path.exists(ifn):
            return ifn
    
    return None
        
    

def build_script(gse):
    fqlist = f'/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/fqlist/{gse}.fq_list.txt'
    
    n_sam = 2  # each project, select n_sam  fastq file for findadapt
    
    fls = []
    cmd_list = []
    
    difficult_gse = {'GSE103493',
        'GSE112289',
        'GSE113411',
        'GSE122488',
        'GSE129255',
        'GSE130654',
        'GSE136321',
        'GSE137617',
        'GSE138983',
        'GSE144781',
        'GSE148576',
        'GSE155281',
        'GSE157528',
        'GSE157916',
        'GSE163507',
        'GSE167863',
        'GSE53451',
        'GSE93020',
        'GSE94721'}
    
    with open(fqlist) as f:
        # GSE146880	/gpfs23/scratch/h_vangard_1/chenh19/exRNA/test/download/SRR11296034/SRR11296034.fastq.gz
        n = 0
        
        
        for i in f:
            if not i.strip():
                continue

            fn_fq_raw = i.strip().split('\t')[1]
            fn_fq = get_fq(fn_fq_raw)
            
            if not fn_fq:
                print(f'fq not found: {fn_fq_raw}')
                continue
            n += 1
            if gse not in difficult_gse:
                if n > n_sam:
                    break
            fls.append(fn_fq)

    if len(fls) < 2:
        print(f'ERROR, fastq file number != 2: {gse}')

    for fq in fls:
        srx_lb = fq.rsplit('/', 2)[-2]
        lb = f'findadapt.{gse}.{srx_lb}'
        cmd = f'fno=/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/realdata/bench_res/{lb}.txt;if [[ ! -s "$fno" ]];then findadapt -fq {fq} -bench &> $fno;else echo "done: {lb}";fi'

        cmd_list.append(cmd)
    
    return cmd_list[-10:]


def main():
    gselist = get_gse_list()
    
    cmds = []
    for gse in gselist:
        tmp = build_script(gse)
        cmds.extend(tmp)
    
    print(len(cmds))
    with open(f'/gpfs23/scratch/h_vangard_1/chenh19/findadapt_bench/realdata/findadapt_cmds.sh', 'w') as o:
        print('\n'.join(cmds), file=o)

main()
