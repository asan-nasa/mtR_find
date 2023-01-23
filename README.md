# mtR_find
## Purpose:

mtR_find is a tool for identification and annotation of sequences mapping to mitochondrial genomes. To run the script on adapter-trimmed FASTQ files, the following dependencies are required: PYTHON pandas module, cutadapt, and bowtie.

## Usage:

mtR_find.py <species_name> <RNA_type> [--FASTA path/to/mitochondrial genome.fa file] [-GTF path/to/gtf file] [--graphical_output yes/no][--output_path path/to/folder] [--input_path path/to/folder] [--files list of files] [--parallel YES/NO] 

The species name and RNA type are required arguments. The others listed above are optional arguments

Valid species name values:<br />
(1) dre - for zebrafish <br />
(2) hsa - for humans <br />
(3) mmu - for mouse <br />
(4) dme - for Drosophila <br />
<br />
Valid RNA type: <br />
(1) sRNA - for mitochondrial sRNA <br />
(2) lncRNA - for mitochondrial long non-coding RNA <br />
<br />
Optional Arguments:<br />
<br />
## Usage examples:
For Zebrafish as species name and if the current directory has all FASTQ files
python mtR_find.py dre
For Zebrafish as species name and to specify filenames explicitly
python mtR_find.py dre --files filename1.fastq filename2.fastq
For Zebrafish as species name and to specify path to folder containing FASTQ files
python mtR_find.py dre --path path/to/folder
For Zebrafish as species name and to specify no graphical output
python mtR_find.py dre --graphical_output no
For Zebrafish as species name and to suspend multiprocessing 

