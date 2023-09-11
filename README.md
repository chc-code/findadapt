# Motivation   
Adapter trimming is an essential step for analyzing small RNA sequencing data where read length is greater than that of target DNA ranging from 18 to 30 bp.   
   
Most adapter trimming tools require adapter information as the input. Adapter information, however, are hard to access, specified incorrectly, or not provided in publicly available datasets, hampering their reproducibility and reusability. Identification of adapter patterns from raw reads is labor-intensive and error-prone. Moreover, the use of randomized adapters to reduce ligation biases in the library preparation makes adapter detection even more challenging. Here, we present FindAdapt, a python package for fast and accurate detection of adapter patterns without any prior information.   
   
# About FindAdapt   
FindAdapt is a python package for identifying and trimming (if specified in the argument) adapter sequences for small RNA sequencing data.   
FindAdapt can also identify the random sequences (such library prepared by NextFlex kit or template switching kits such as SMARTer / CATS kit).   
FindAdapt do not depend on prior knowledge about the adapter, so it can also correctly identify the customized adapter sequences.   
   
We also tested the package on fastq files from 3148 biosamples which covers different library prepartion kits, adapter pattern and reads structure, FindAdapt can correctly identify the adapter sequences.   
   
FindAdapt is written in Python3 and is a stand-alone script without any package requirements (If need to trim fastq file after adapter detection, you need to install cutadapt package first).   
   
# Usage   
   
The usage can be accessed by   
`findadapt -h`   
If your input data are human small RNA seq, you just need to specify the parameters in the `Main arguments` section below   

## Main arguments    
You can only select one option (either `-fq` or `-prj`) as the input

- `fn_fq_file`   Optional positional argument,  the path for single fastq file   
- `-fn_fq_list / -list / -l file_list`   text file, contains the fastq list. column1 = study ID, column2 = fastq file path.  -fq and -prj are mutually exlusive   

## Optional arguments

###  customize reference   
By default, the reference is most abundant 100 miRNA in human. you can also useyour own sequence list if the your input doesn't match the default reference setting, available options are   

- `-fn_refseq filename`  User specified sequence known to express in the dataset, can be fasta format or one sequence per line.   
- `-organism / -org str` organism name  default is human. valid = human, mouse, fruitfly, worm (c. elegans), arabidopsis, rice. Alternatively, you can use the miRBase prefix, such hsa, mmu, dme, cel, ath, osa; or, if -fn_refseq is specified, you can specify as other
- `-list_org`  list the supported organism and exit


### Output options   
- `-o prefix`,str, optional, the prefix for the output results, if not specified, will infer from the input file
- `-quiet / -q` , toggle, suppress the pyahocorasick not installed warning
- `-cut / -cutadapt/ -trim`  flag,  run the cutadapt process, need the cutadapt already installed and available in PATH   
- `-pw_cutadapt str`  specify the cutadapt path, default is search from PATH   
- `-v / -verbose` flag, display the full logging information in the terminal   
   
   
### Run control

- `-expected_adapter_len int`  expected adapter length to discover, by default, 12 nt   
- `-max_random_linker int`   max allowed random seqence length, default = 6
- `-nreads int`  max reads number used to find adapter, default is 1 million, if use all reads, set as -1
- `-nsam int`  for studies with multiple samples/fastq files, by default, will use first 5 files to infer the adapter pattern for this study. You can change the number by this parameter. If you need to use all samples, set as -1   
- `-thres_multiplier float`    the multiplier threshold for accepting an extra random base from a lower phase. e.g for 2 possible (3p-random_sequence_lengh, adapter_seq) : (0, seq1), (1, seq2), the reads supporting each combination is n1 and n2 respectively,  (1, seq2) will only be accepted if n2 > n1 * thres_multiplier. this argument is to avoid the including 

- `-min_reads int`  minimum matched reads number for infering per fastq file, default=30, if lower than this value, the adapter inferring step will be skip, you may need to check the reference settings.   
- `-threads / -cpu int` the threads to use, by default = 5. Usually, the performance won't improve greatly when more than 5 threads are used (refer to the manuscript for detail)   
- `-enough_reads int` enough matched reads number for infering per fastq file, after reaching, will stop reading the raw fastq file, default=1000   
- `-f`  `-force`  flag, force rerun the analysis, ignoring the exisiting parsed reads,  can be useful when you used new reference.   


# Examples   
   
Here we provided several fastq files from 3 studies   
1. GSE106303, the adapter sequence is not specified in the GEO database or article   
2. GSE122068, NextFLEX library preparation kit, the reads will have 4N random sequence at both 5' and 3' end of the insert   
3. GSE137617, SMARTer library preparation kit, the 3' adapter sequence of the reads is polyA, and have multiple (usually 3 nt) random sequence at the 5' end   
   
you can run get the adapter by   
`./findadapt -fq <fn_fq>`   
   
 for example, GSE122068.nextflex.SRR8144939.truncated.fastq.gz   
 `./findadapt -fq ./demo/GSE122068.nextflex.SRR8144939.truncated.fastq.gz`   
