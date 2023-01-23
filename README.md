# mtR_find
## Purpose:

mtR_find is a tool for identification and annotation of sequences mapping to mitochondrial genomes. To run the script on adapter-trimmed FASTQ files, the following dependencies are required: PYTHON pandas module, cutadapt, and bowtie.

## Usage:

mtR_find.py <species_name> <RNA_type> [--FASTA path/to/mitochondrial genome.fa file] [-GTF path/to/gtf file] [--graphical_output yes/no][--output_path path/to/folder] [--input_path path/to/folder] [--files list of files] [--parallel YES/NO] 

Required arguments: The species name and RNA type are required arguments. The others listed above are optional arguments

Valid species name values:<br />
(1) dre - for zebrafish <br />
(2) hsa - for humans <br />
(3) mmu - for mouse <br />
(4) dme - for drosophila <br />
(5) xen - for xenophus <br />
(6) gal - for chicken <br />
(7) rno - for rat <br />
(7) non_model - for species not listed above <br />

<br />
Valid RNA type: <br />
(1) sRNA - for mitochondrial sRNA <br />
(2) lncRNA - for mitochondrial long non-coding RNA <br />
<br />
Optional Arguments:<br />
<br />
--parallel:  default value is "NO". If users want to suspend multuprocessing, they have to specify "YES".<br />
--FASTA and -GTF: by default if users specify anyone of the 6 species code (listed above), the script would download the FASTA and GTF file automatically. In case if users want to analyze mtsRNAs/mtlncRNAs in any other species, they would have to manually download the mitochondrial genome FASTA file
<br />
--input_path: defalut value is current working directory. But users can specify a input path<br />
--graphical_output: No. If graphical output of the basic plots has to be generated, the user has to specify
“yes” under graphical_output <br />
--files: defalut value is "None". If the files are in different locations, the absolute path of the files can be
specified. Note: if –input_path is specified, --files cannot be specified. If –files
argument is specified, files in the current working directory will not be analyzed,
even though the output directory will be the current working directory – unless a
different output path is specified using –output_path argument. <br />

## Usage examples:
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
python mtR_find.py dre --parallel NO  <br />
<br />
To analyze mtsRNA with zebrafish as species name and if the current directory has all FASTQ files <br />
python mtR_find.py dre lncRNA  <br />

