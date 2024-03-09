# mtR_find
## Purpose:

mtR_find is a tool for identification and annotation of sequences mapping to mitochondrial genomes. To run the script on adapter-trimmed FASTQ files, the following dependencies are required: pandas (version 0.21.0 and above), multiprocessing, bowtie (version 1.1.2 and above) and samtools (version 1.9 and above). If users want to output basic plots, then matplotlib PYTHON module is also needed (optional - version 2.0.2 and above). 

## Usage:

mtR_find.py <species_name> <RNA_type> [--FASTA path/to/mitochondrial genome.fa file] [-GTF path/to/gtf file] [--graphical_output yes/no][--output_path path/to/folder] [--input_path path/to/folder] [--files list of files] [--parallel YES/NO] 

Required arguments: The species name and RNA type are required arguments. The others listed above are optional arguments

Valid species name values:<br />
(1) dre - for zebrafish <br />
(2) hsa - for humans <br />
(3) mmu - for mouse <br />
(4) xen - for xenophus <br />
(5) gal - for chicken <br />
(6) rno - for rat <br />
(7) non_model - for species not listed above <br />

<br />
Valid RNA type: <br />
(1) sRNA - for mitochondrial sRNA <br />
(2) lncRNA - for mitochondrial long non-coding RNA <br />
<br />
Optional Arguments:<br />
<br />
--parallel:  default value is "NO". If users want to suspend multuprocessing, they have to specify "YES".<br />
--FASTA and -GTF: by default if users specify anyone of the 7 species code (listed above), the script would download the FASTA and GTF file automatically. In case if users want to analyze mtsRNAs/mtlncRNAs in any other species, they would have to manually download the mitochondrial genome FASTA file
<br />
--input_path: defalut value is current working directory. Users can specify a input path <br />
--output_path: defalut value is current working directory. Users can specify a output path <br />
--cutoff: default value = 200. cutoff corresponds to the threshold value of the total ncRNA count of individual ncRNAs from all the libraries combined together. For example, a cutoff value of 200 would discard ncRNA sequences with total count (from all libraries) less than 200.<br />
--filter, default = 50, a length filter applicable only for mt-lncRNAs, if users want to study only lncRNAs greater than 200, they can specify "--filter 200" in the command line<br />
<br />
--files: defalut value is "None". If the files are in different locations, the absolute path of the files can be
specified. Note: if –input_path is specified, --files cannot be specified. If –files
argument is specified, files in the current working directory will not be analyzed,
even though the output directory will be the current working directory – unless a
different output path is specified using –output_path argument. <br />
--graphical_output: default values is "no". If graphical output of the basic plots has to be generated, the user has to specify “yes” under graphical_output. If "YES" is specified, then it is also mandatory to specify a metadata file. This metadata file should contain the filenames as the first column and the condition/factor as the second column <br />

## Usage examples:
To analyze mtsRNA with human as species name and if the current directory has all FASTQ files <br />
python mtR_find.py hsa sRNA <br />
To analyze mtsRNA with zebrafish as species name and if the current directory has all FASTQ files <br />
python mtR_find.py dre sRNA <br />
<br />
To analyze mtsRNA with zebrafish as species name and to specify filenames explicitly <br />
python mtR_find.py dre --files filename1.fastq filename2.fastq <br />
<br />
To analyze mtsRNA with zebrafish as species name and to specify path to folder containing FASTQ files <br />
python mtR_find.py dre --path path/to/folder <br />
<br />
To analyze mtsRNA with zebrafish as species name and to specify no graphical output <br />
python mtR_find.py dre --graphical_output no <br />
<br />
To analyze mtsRNA with zebrafish as species name and to suspend multiprocessing <br />
python mtR_find.py dre --suspend yes  <br />
<br />
To analyze mtsRNA with zebrafish as species name and if the current directory has all FASTQ files <br />
python mtR_find.py dre lncRNA  <br />

Note: bowtie and samtools should be added to PATH. More information can be found here:<br />

https://phoenixnap.com/kb/linux-add-to-path <br />
<br />
https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path <br />
<br />
https://opensource.com/article/17/6/set-path-linux <br />

# Test datasets
## Adapter trimmed FASTQ files:

Adpater trimming was performed on test datsets using Cutadapt version 1.5. If users are trimming raw files downloaded from SRA using other adpater trimming tools (other than Cutadapt) or differnt version of Cutadapt, there could be variations to the input data and hence would affect reproducibility of the results. To ensure reprdocubility, we have uploaded the adapter trimmed FASTQ files for test datasets and can be downloaded in the links provided below. Interested users who wanted to test mt_find on test data before testing on their own data, can download the adpater trimmed FASTQ files from the link below.  

## Dataset-1:
https://filesender.sikt.no/?s=download&token=dcb87ce5-7224-4900-bd3b-5ec1a4295d09
## Dataset-2:
https://filesender.sikt.no/?s=download&token=21680929-cd29-47d7-839c-91c5cf992052
## Dataset-3:
https://filesender.sikt.no/?s=download&token=7547e686-11a8-4c62-a19b-67c888a6b904

## Results:

The results from running mtR_find for each test dataset which inlcudes the output from mtR_find and nohup commonad line output can be found in the test folder, under sub-folders dataset-1, dataset-2 and dataset-3 respectively
