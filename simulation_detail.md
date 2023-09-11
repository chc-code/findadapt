# 1. original data used for simulation

The reads for simulation come from /scratch/h_vangard_1/chenh19/findadapt_bench/simudata/SRR6502962.trimmed.fastq.gz

kit = NEBNext

The above file is got from running

cutadapt -j4 -q 20 -a AGATCGGAAGAG -O 1 --trim-n -m 16 -o SRR6502962.trimmed.fastq.gz /gpfs23/scratch/h_vangard_1/chenh19/exRNA/test/download/SRR6502962/SRR6502962.fastq.gz

below is the log for the cutadapt run
```bash
         Total reads processed:                 963,413
         Reads with adapters:                   958,326 (99.5%)
         == Read fate breakdown ==
         Reads that were too short:              36,330 (3.8%)
         Reads written (passing filters):       927,083 (96.2%)
         Total basepairs processed:    48,170,650 bp
         Quality-trimmed:                  57,405 bp (0.1%)
         Total written (filtered):     25,703,471 bp (53.4%)
```

# 2. parsing the reads

Parse the above SRR6502962.trimmed.fastq.gz (total reads = 927083 ) into a dict,
Key = index (from 0 to 927083)

# 3. build simulation dataset

## 3.1 select the 3’ adapter sequence

If the setting is  reads_struc type 1-4  (ordinary 3’ adapter),  will randomly pick an adapter from the following adapter pool,

```
    'AGATCGGAAGAG',  # NEBNext
    'AACTGTAGGCAC', # kiaseq
    'TGGAATTCTCGG', # nextflex, truseq
```

If the setting is reads_struc type 5-8, will use 12nt poly(A)  as the 3’ adapter
If the settings is no adapter, the 3’ adapter sequence is empty string
This adapter sequence will be used throughout the whole dataset.

## 3.2 sample the trimmed reads

selected_reads_index  =  sample 10k items from the key of the reads dict
## 3.3 build the new reads

Iterate over above selected_reads_index
in each loop: `original sequence = reads_dict[idx]`

If add_5p_randomer or add_3p_randomer is True, will

1. random choose a number from 1-5 (len_randomer)
2. random choose from `[A, T, C, G]`  len_randomer times to build the random-mer sequence

the quality score will be filled with “F” for the added sequence (random-mer or adapter sequence)

new read sequence =
5p randomer sequence + original sequence + 3’ randomer + adapter sequence
