# mtR_find
## Purpose:

mtR_find is a tool for identification and annotation of sequences mapping to mitochondrial genomes. To run the script on adapter-trimmed FASTQ files, the following dependencies are required: PYTHON pandas module, cutadapt, and bowtie.

## Usage:

mtR_find.py <species_name> <RNA type> [--FASTA path/to/mitochondrial genome.fa file] [-GTF path/to/gtf file] [--graphical_output yes/no][--output_path path/to/folder] [--input_path path/to/folder] [--files list of files] [--parallel YES/NO] 
