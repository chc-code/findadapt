# Motivation   
Adapter trimming is the first step for analyzing small RNA sequencing data where reads are longer than target RNAs with lengths ranging from 18 to 30 bp. There is a lack of tools for accurately identifying adapters from raw reads. Moreover, the use of randomized adapters to reduce ligation biases in small RNA-seq library preparation makes adapter detection even more challenging. 
   
# About FindAdapt   
FindAdapt is a Python package for identifying adapter patterns for small RNA sequencing data without dependency on prior information.   

# Installation
FindAdapt is a stand-alone Python package (python >=3.6).

## Download and uncompress
```bash
wget https://github.com/chc-code/findadapt/archive/refs/heads/master.zip
unzip master.zip  # the output folder will be findadapt-master

# use FindAdapt
cd findadapt-master
./findadapt  -h
```

To achieve the best performance, the installation of pyahocorasick is highly recommended. 
```
# install pyahocorasick
pip install pyahocorasick
```

## docker image
A docker image is also available at https://hub.docker.com/r/chccode/findadapt
(`pyahocorasick` is contained, and the findadapt script is set as the entrypoint.)

```
docker pull chccode/findadapt
docker run  chccode/findadapt reads.fastq.gz  
```

 
# SYNOPSIS   

      findadapt reads.fastq.gz
      findadapt reads.fastq.gz -organism mouse
      findadapt -list_org
      
      # identify adapter sequences for all samples of a study from mice 
      findadapt -fn_fq_list fq_list.txt -organism mouse

      # trim the identified adapters using cutadapt package 
      findadapt reads.fastq.gz -pw_cutadapt path/to/cutadapt -cut  


# COMMANDS AND OPTIONS  

## File input
You can only select one option (either `-fq` or `-prj`) as the input

- `fn_fq_file`   Optional positional argument,  the path for single fastq file   
- `-fn_fq_list / -list / -l file_list`   a tab-delimited file, containing the list of fastq files. column1 = study ID, column2 = path of the fastq file. 

## Reference sequences
either a list of sequences (fasta format or one sequence per line) by '-fn_refseq' or organism name by `-organism`
- `-fn_refseq filename`   a list of sequences in fasta format or one sequence per line.   
- `-organism / -org str`   organism name (for example, human, mouse, fruitfly, worm, arabidopsis, rice);  default: human.  Alternatively, the miRBase prefix, such hsa, mmu, dme, cel, ath, osa; or, if -fn_refseq is specified, you can specify as other
- `-list_org`  list the supported organisms 


## Output options   
- `-o prefix`,str, optional, the prefix for the output results, if not specified, will infer from the input file
- `-quiet / -q` , toggle, suppress the warning message if pyahocorasick not installed 
- `-cut / -cutadapt/ -trim`  flag,  run the cutadapt process; require the cutadapt already installed and available in PATH   
- `-pw_cutadapt str`  the path of cutadapt, the default is from PATH   
- `-v / -verbose` flag, display the log information in the terminal   
   
   
## Other Options
- `-expected_adapter_len int`  the length of adapter sequence, default = 12 bp   
- `-max_random_linker int`   the maximum length of random-mer, default = 8 bp
- `-nreads int`  the maximum number of reads used to find adapter, default: 1 million, if use all reads, set as -1
- `-nsam int`  the number of samples foradapter identification in a file list,  default: all samples
- `-thres_multiplier float`    the threshold of the ratio between the count of the child and the count of the parent, default=1.2; if >=1.2, save the child record; otherwise, save the parent record
- `-min_reads int`  the minimum number of matched reads for adapter identification, default=30. if lower than this value, the adapter identification will fail and users may need to check the reference settings.   
- `-threads / -cpu int` the number of threads, default = 5. 
- `-enough_reads int` the number of matched reads for adapter identification, default=1000   
- `-f`  `-force`  flag, force rerun the analysis, ignoring the exisiting parsed reads,  can be useful when use a new reference.   


# Examples   
   
We provided several fastq files from three studies   
1. GSE106303, the adapter sequence is not specified in the GEO database or the literature  
2. GSE122068, generated by NextFLEX library preparation kit where reads  have 4N random sequence at both the 5' and 3' ends 
3. GSE137617, generated by SMARTer library preparation kit where multiple (usually 3 nt) random bps at the 5' end and polyA as the 3' adapter sequence
   
To identify adapter sequences  
`./findadapt <fn_fq>`   
   
 for example, GSE122068.nextflex.SRR8144939.truncated.fastq.gz   
 `./findadapt ./demo/GSE122068.nextflex.SRR8144939.truncated.fastq.gz`   

 # Output Format
 ## log information
 ```
 2023-09-08 08:13:02  INFO   <module>              line: 1683   1/1: single - using 1/ 1 fq files
2023-09-08 08:13:02  INFO   get_adapter_per_prj   line: 1076   	processing GSE122068.nextflex.SRR8144939.truncated.fastq.gz
2023-09-08 08:13:02  INFO   get_parsed_reads      line: 834    matched reads found: 1177
2023-09-08 08:13:02  INFO   export_data           line: 1229   	most possible kit = NEXTflex
2023-09-08 08:13:02  INFO   export_data           line: 1289   result per-prj = GSE122068.nextflex.SRR8144939.truncated.adapter.txt
2023-09-08 08:13:02  INFO   export_data           line: 1290   result per-fq = GSE122068.nextflex.SRR8144939.truncated.per_fq.adapter.txt
```
## .adapter.txt
The output contains the following columns
**Prj**: The output prefix, if the input is a single fastq file rather than a fastq file list (-fn_fq_list), it will be "single"
**total_reads**: Total matched reads used for adapter identification
**3p_seq**:  The sequence of 3' adapter
**3p_phase**: the random sequence length before 3' adapter
**3p_count / 3p_reatio**: The number and ratio of reads supporting this 3' adapter sequence and random sequence length
**5p_phase**:  the random sequence length before the insert
**5p_count / 5p_ratio**:  The number and ratio of reads supporting this 5' random sequence length 
**err**:  the error information if fail to get the adapter sequence, 

You can remove the adapter using the identified pattern by specifying `-cut` 
Or use the output to build your own cutadapt command.

```
# if 3p_seq is empty and 5p_phase > 0:
{pw_cutadapt} -u {5p_phase} -m 15 -j 8  --trim-n {fn_fq} -o {fn_out}

# elif 3p seq is not empty and 5p_phase = 3p_phase = 0
{pw_cutadapt} -a {seq_3p} -m 15 -j 8  --trim-n  {fn_fq} -o {fn_out}

# if 3p_phase > 0 and 5p_phase == 0
{pw_cutadapt} -a {seq_3p} -j 8 --trim-n  {fn_fq} |cutadapt -u -{3p_phase} -m 15 -o {fn_out}

# if 3p_phase = 0 and 5p_phase > 0
{pw_cutadapt} -a {seq_3p} -j 8 --trim-n  {fn_fq} |cutadapt -u {5p_phase} -m 15 -o {fn_out}

# if 3p_phase > 0 and 5p_phase > 0
{pw_cutadapt} -a {seq_3p} -j 8 --trim-n  {fn_fq} |cutadapt -u -{3p_phase} -u {5p_phase} -m 15 -o {fn_out}
```


| prj    | total_reads | 3p_seq       | 3p_phase | 3p_count | 3p_ratio | 5p_phase | 5p_count | 5p_ratio | err |
|--------|-------------|--------------|----------|----------|----------|----------|----------|----------|-----|
| single | 1177        | TGGAATTCTCGG | 4        | 1021     | 0.8667   | 4        | 1143     | 0.9711   |


## .per_fq.adapter.txt
The detail adaoter information of every input fastq file

| prj    | fastq                                   | total_reads | side | sn | seq          | phase | count | ratio  |
|--------|-----------------------------------------|-------------|------|----|--------------|-------|-------|--------|
| single | GSE122068.nextflex.SRR8144939.truncated | 1177        | 3p   | 1  | TGGAATTCTCGG | 4     | 1021  | 0.8675 |
| single | GSE122068.nextflex.SRR8144939.truncated | 1177        | 3p   | 2  | CTGGAATTCTCG | 3     | 633   | 0.5378 |
| single | GSE122068.nextflex.SRR8144939.truncated | 1177        | 5p   | 1  |              | 4     | 1143  | 0.9711 |

