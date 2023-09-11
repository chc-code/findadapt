"""
build the simulation fastq file for performance test

acutal fastq file used is
/scratch/h_vangard_1/chenh19/exRNA/rerun/download/SRR13873945/SRR13873945.fastq.gz
which is 4N+4N pattern
"""
import re, sys, os, gzip

import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('fn', help="""input fastq file""")
ps.add_argument('-force', '-f', action='store_true', help="""ignore the prev data""")
args = ps.parse_args()

sys.path.append('/home/chenh19/jb/query')

from findadapt import get_most_expressed_smallRNA

def getlogger(prefix=None):
    try:
        # return the logger defined in this script
        return globals()['logger']
    except:
        pass

    import sys, logging
    logger_name = __file__
    try:
        logger = logging.getLogger(logger_name)
    except:
        logger = logging.getLogger('terminal')
    fmt = logging.Formatter('%(asctime)s  %(levelname)-6s %(funcName)-25s  line: %(lineno)-5s  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logger.setLevel('DEBUG')

    handler_names = {_.name for _ in logger.handlers}
    if 'console' not in handler_names:
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(fmt)
        console.setLevel('INFO')
        console.name = 'console'
        logger.addHandler(console)

    if prefix and 'file' not in handler_names:
        fn_log = f'{prefix}.log'
        fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
        fh_file.setLevel('DEBUG')
        fh_file.setFormatter(fmt)
        fh_file.name = 'file'
        logger.addHandler(fh_file)

    return logger

logger = getlogger()

def prepare_data(fn, force=False):

    fn_matched = 'matched_reads.fastq'
    fn_unmatched = 'unmatched_reads.fastq'
    if not force and os.path.exists(fn_matched) and os.path.exists(fn_unmatched):
        logger.info(f'raw data already done')
        return 0

    fh_matched = open(fn_matched, 'w')
    fh_unmatched = open(fn_unmatched, 'w')



    matched = {'max': 1_000_000, 'ct': 0, 'fh': fh_matched, 'lb': 'matched', 'done': 0}
    unmatched = {'max': 2_000_000, 'ct': 0, 'fh': fh_unmatched, 'lb': 'unmatched', 'done': 0}

    res = {0: unmatched, 1: matched}

    mir_list = get_most_expressed_smallRNA('human')
    mir_list = [_[0] for _ in mir_list]

    def check_match(seq):
        for mir_seq in mir_list:
            if mir_seq in seq:
                return 1
        return 0

    logger.info('now reading the fq file')
    with gzip.open(fn, 'rt') as f:
        n = 0
        while True:
            if res[0]['done'] + res[1]['done'] > 1:
                logger.info(f'enough reads got, n = {n}')
                break
            header = f.readline()
            if not header:
                logger.warning(f'file end reached')
                break
            if header[0] != '@':
                continue

            n += 1
            seq = f.readline()
            dummy = f.readline()
            qual = f.readline()

            if n % 100000 == 0:
                logger.info(f'processing {n/10000:.0f}w reads')
            flag = check_match(seq)
            ires = res[flag]
            if not ires['done'] and ires['ct'] < ires['max']:
                fh = ires['fh']
                ires['ct'] += 1
                fh.write(header)
                fh.write(seq)
                fh.write(dummy)
                fh.write(qual)
            else:
                ires['done'] = 1

    for v in res.values():
        lb = v['lb']
        logger.info(f'{lb}:  ct = {v["ct"]}')
        v['fh'].close()




def build_data():
    fn_matched = 'matched_reads.fastq'
    fn_unmatched = 'unmatched_reads.fastq'

    fh_matched = open(fn_matched, 'r')
    fh_unmatched = open(fn_unmatched, 'r')


    def write_data(flag, fho):
        fh = fh_matched if flag else fh_unmatched
        for _ in range(4):
            fho.write(fh.readline())

    import random
    n_total_exp = 1_500_000


    for ratio in [0.9]:
    # for ratio in [0.5, 0.1, 0.01, 0.001]:
        logger.info(f'building data for ratio {ratio}')
        n_matched_written = 0
        n_unmatched_written = 0
        n_total = 0
        fn = f'benchmark_test.miRNA_ratio_{ratio}.fastq'
        fho = open(fn, 'w')
        fh_matched.seek(0)
        fh_unmatched.seek(0)

        while True:
            _ = random.random()

            if n_total > n_total_exp:
                logger.info(f'matched written = {n_matched_written}, n_unmatched_written = {n_unmatched_written}, , total = {n_total}, actual ratio = {n_matched_written/n_total:g}')
                break

            if _ < ratio:
                n_matched_written += 1
                n_total += 1
                write_data(1, fho)
            else:
                n_total += 1
                n_unmatched_written += 1
                write_data(0, fho)

        fho.close()

    fh_matched.close()
    fh_unmatched.close()







if __name__ == "__main__":
    fn = args.fn
    force = args.force
    prepare_data(fn, force=force)
    build_data()
