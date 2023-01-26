#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
import pandas as pd
import os
import sys
from collections import Counter
import operator
from itertools import takewhile
import multiprocessing
from functools import partial
import argparse
import itertools

a = int(pd.__version__.split(".")[0])
b = int(pd.__version__.split(".")[1])
if a >= 1:
   print("Pandas version is " + str(pd.__version__) + " check Okay!")
elif b>= 21:
   print("Pandas version is " + str(pd.__version__) + " check Okay!")
else:
   sys.exit('\nERROR: mtR_find requires Pandas version 0.21.0 or higher. Please update Pandas\n')

docstring= """

USAGE:
python mtR_find.py <species_name> <sRNA/lncRNA>

Arguments:
Valid species name
(1) dre - for zebrafish
(2) hsa - for humans
(3) mmu - for mouse
(4) dme - for Drosophila
(5) xen - for xenophus
(6) gal - for chicken
(7) rno - for rat
(7) non_model - for species not listed above

Valid RNA type:
(1) sRNA
(2) lncRNA
Use --help for more info

DESCRIPTION
Creates a read count of mitochondrial derived sequences and also outputs the annotation information

"""

parser = argparse.ArgumentParser(usage = "\n" +"python %(prog)s species_name [--non_model NON_MODEL] [--FASTA FASTA] [-GTF GTF] [--graphical_output GRAPHICAL_OUTPUT][--output_path path/to/folder] [--input_path path/to/folder] [--files list of files] \n" +"\n" + "Usage examples:\nFor Zebrafish as species name and if the current directory has all FASTQ files\npython mtR_find.py dre\nFor Zebrafish as species name and to specify filenames explicitly\npython mtR_find.py dre --files filename1.fastq filename2.fastq\nFor Zebrafish as species name and to specify path to folder containing FASTQ files\npython mtR_find.py dre --path path/to/folder\nFor Zebrafish as species name and to specify no graphical output\npython mtR_find.py dre --graphical_output no\n" + "\n" + "Description:\nCreate read count file with sequences mapping to mitochondrial genome with their annotation information\n", add_help=False)
required = parser.add_argument_group('required arguments')
required.add_argument("species_name", help = """enter species name:Accepted species name arguments: (1) species name = dre #(for zebrafish) (2) species name = hsa #(for humans)(3) species name = mmu #(for mouse) (4) species name = dme #(for Drosophila) (5) species name = xen #(for xenophus) (6) species name = gal #(for chicken)""")
required.add_argument("RNA", help = """enter sRNA for small non-coding RNA and lncRNA fro long non-coding RNA""")
optional = parser.add_argument_group('optional arguments')
optional.add_argument("--non_model", help="""for non-model organisms or when the mitochondrial genome/annotation file is already downlaoded locally""")
optional.add_argument("--FASTA", help = "specify the absolute path to the mitochondrial genome FASTA file")
optional.add_argument("--GTF", help = "specify the absolute path to the mitochondrial genome GTF file")
optional.add_argument("--cutoff", default = 200, help = "filter sRNAs that have a combined total count below a specified value in all libraries")  

optional.add_argument('--input_path', default = os.getcwd(), help= 'paste path to FASTQ files')
optional.add_argument('--output_path', default = os.getcwd(), help= 'paste path to store output files')
optional.add_argument('--suspend', default = "no", help= 'option to suspend multiprocessing')
optional.add_argument('--files', nargs='*', help= 'enter FASTQ files seperated by space(enter the absolute path and not the relative path)')

optional.add_argument("--filter", default = 51, help = "applicable only for mtlncRNA, if this option is specified then lncRNAs above --filter value are reported")
optional.add_argument("--graphical_output", default = "no", help="""If you want to specify graphical output specify -grapical_output yes in the command line followed by species name""")
optional.add_argument("--metadata", help= 'enter PATH to metadata file')
#color optional only if the number of elements in the condition/factor exceeds 16
optional.add_argument("--color", nargs='*', help= 'enter color names entered by space')
optional.add_argument("-h", "--help", action='help', help='print help message')
args = parser.parse_args()

if args.graphical_output == "yes":
 try:
   import matplotlib
 except ImportError:
  print("Install matplotlib")
 else:
  a = matplotlib.__version__
  if float(a[:3]) >= 2:
      print("Matplotlib version check OK......")
  else:
      sys.exit('\nERROR: mtR_find requires Matplotlib version 2.0.2 or higher\n')
 if args.metadata == None:
    sys.exit('\nERROR: Metadata option is mandatory if graphical_output == "yes"\n')

if (args.suspend != "yes") and (args.suspend != "no"):
   sys.exit('\nERROR: --suspend ca take only either one of the arguments: "yes" or "no"\n')

mt_list=['tRNA-Phe', 'mtSSU rRNA', 'tRNA-Val', 'mtLSU rRNA', 'tRNA-Leu', 'ND1', 'tRNA-Ile', 'tRNA-Gln', 'tRNA-Met', 'ND2', 'tRNA-Trp', 'tRNA-Ala', 'tRNA-Asn', 'tRNA-Cys', 'tRNA-Tyr', 'COI', 'tRNA-Ser1', 'tRNA-Asp', 'COII', 'tRNA-Lys', 'ATP8', 'ATP6', 'COIII', 'tRNA-Gly', 'ND3', 'tRNA-Arg', 'ND4L', 'ND4', 'tRNA-His', 'tRNA-Ser2', 'tRNA-Leu2', 'ND5', 'ND6', 'tRNA-Glu', 'CytB', 'tRNA-Thr', 'tRNA-Pro']

#check if bowtie exist
os.system("bowtie --version > vers.txt")
infile= open("vers.txt", "r")
lines = infile.readlines()
if lines != []:
   ver = lines[0].strip().split(" ")[2]
   print("Bowtie version is  " + str(ver) + " check okay!")
else:
   ver = "no"
if (ver == "no"):
   sys.exit('\nERROR: Bowtie not found. Please add the PATH to bowtie $PATH\n%s')
else:
   bowtie = "bowtie"


def worker(f):
    infile= open(f)
    filename = str(f.split(".")[0].split("_")[0])
    fastq_lst = infile.readlines()[1::4]
    d = len(fastq_lst)
    fastq_lst = [line.strip() for line in fastq_lst]
    b = Counter(fastq_lst)
    return d,filename,b


#defining how to extract MT cordinates from gtf file
def extract_MT(infile):
    print("Extracting mitochondrial gene annotation info.....")
    lines = infile.readlines()[5:]
    annotation=[]
    i = 0
    for line in lines:
     line= line.strip().split("\t")
     if (args.species_name != "xen"):
      if (line[0] == 'MT' and line[2]== 'gene'):
        feature=line[8].split(";")
        biotype=feature[4].split(" ")[::-1][0].replace('"',"")
        annotation.append((line[3],line[4],line[6],mt_list[i],biotype))
        i = i + 1
     else:
      if (line[0] == 'MT') and (line[2]== 'gene' or line[2]== 'tRNA' or line[2]== 'rRNA'):
        feature=line[8].split(";")
        if (line[2] == "gene"):
           biotype=feature[4].split("=")[1]
           annotation.append((line[3],line[4],line[6],feature[-1].split("=")[1],biotype))
        else:
           biotype=feature[-1].split("=")[1]
           annotation.append((line[3],line[4],line[6],feature[-2].split("=")[1],biotype))
        i = i + 1

    return annotation

#define function to annotate fragments
def mt_annotator(line,kind):
 element = []
 line= line.strip().split("\t")
 subs = line[12].split("MD:Z:")[1]
 if (kind != "mt-lncRNA") and len(subs) > 2:
    subs = subs
 elif (kind == "mt-lncRNA") and (len(subs) > 3):
    subs = subs
 else:
    subs = ""
 endposition = (len(line[0])+int(line[3]))-1
 i = 0
 for x in annotation:
  i=i+1
  if(i<37):
   if (int(x[0]) <= int(line[3]) < int(x[1])) and (len(line[0])+int(x[0])<=int(x[1])):
    if line[1] == "0" and x[2] == "+":
      element.extend([line[0],kind,x[4],"H","sense",x[3].split("-")[-1],"within gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "0" and x[2] == "-":
      element.extend([line[0],kind,x[4],"H","anti-sense",x[3].split("-")[-1],"within gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "16" and x[2] == "+":
      element.extend([line[0],kind,x[4],"L","anti-sense",x[3].split("-")[-1],"within gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "16" and x[2] == "-":
      element.extend([line[0],kind,x[4],"L","sense",x[3].split("-")[-1],"within gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
   if int(x[0]) <= (int(line[3])+len(line[0])-1) <= int(x[1]):
    if line[1] == "0" and x[2] == "+":
      element.extend([line[0],kind,x[4],"H","sense",x[3].split("-")[-1],"overlap gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "0" and x[2] == "-":
      element.extend([line[0],kind,x[4],"H","anti-sense",x[3].split("-")[-1],"overlap gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "16" and x[2] == "+":
      element.extend([line[0],kind,x[4],"L","anti-sense",x[3].split("-")[-1],"overlap gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "16" and x[2] == "-":
      element.extend([line[0],kind,x[4],"L","sense",x[3].split("-")[-1],"overlap gene boundary",line[3],endposition,x[0],x[1],subs])
      return element
  else:
    if line[1] != "16":
     element.extend([line[0],kind,"non-coding","H","sense","non-coding","falls in non-coding region",line[3],endposition,"na","na",subs])
     return element
    else:
     element.extend([line[0],kind,"non-coding","L","anti-sense","non-coding","falls in non-coding region",line[3],endposition,"na","na",subs])
     return element
# annotation for long non coding RNA
def mt_annotator2(line,kind):
 element = []
 line= line.strip().split("\t")
 subs = line[12].split("MD:Z:")[1]
 if (kind != "mt-lncRNA") and len(subs) > 2:
    subs = "-" + subs
 elif (kind == "mt-lncRNA") and (len(subs) > 3):
    subs = "-" + subs
 else:
    subs = ""
 endposition = (len(line[0])+int(line[3]))-1
 i = 0
 for x in annotation:
  i=i+1
  if(i<37):
   if (int(x[0]) <= int(line[3]) < int(x[1])):
    if line[1] == "0" and x[2] == "+":
      element.extend([line[0],kind,x[4],"H","sense",x[3].split("-")[-1],"seq start site inside gene",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "0" and x[2] == "-":
      element.extend([line[0],kind,x[4],"H","anti-sense",x[3].split("-")[-1],"seq start site inside gene",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "16" and x[2] == "+":
      element.extend([line[0],kind,x[4],"L","anti-sense",x[3].split("-")[-1],"seq start site inside gene",line[3],endposition,x[0],x[1],subs])
      return element
    if line[1] == "16" and x[2] == "-":
      element.extend([line[0],kind,x[4],"L","sense",x[3].split("-")[-1],"seq start site inside gene",line[3],endposition,x[0],x[1],subs])
      return element

  else:
    if line[1] != "16":
     element.extend([line[0],kind,"non-coding","H","sense","non-coding","seq start site outside gene",line[3],endposition,"na","na",subs])
     return element
    else:
     element.extend([line[0],kind,"non-coding","L","anti-sense","non-coding","seq start site outside gene",line[3],endposition,"na","na",subs])
     return element


#download mitochondrial genome of species and create bowtie index and aslo downlaod gtf file and extract MT coordinates

if (args.species_name == "dre") and (args.non_model==None):
   if not os.path.exists("zebrafish"):
    os.makedirs("zebrafish")
   species = "zebrafish"
   os.chdir("zebrafish")
   if not os.path.exists("MT_genome"):
    os.makedirs("MT_genome")
   os.chdir("MT_genome")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.MT.fa.gz"
   print("Downloading Zebrafish Mitochondrial genome....................")
   print(command)
   os.system(command)
   os.system("gunzip -f *.gz")
   os.chdir("../")
   print("Creating bowtie index for Zebrafish mitochondrial genome")
   if not os.path.exists("bowtie-index"):
    os.makedirs("bowtie-index")
   index = "bowtie-index/zebra_MT_index"
   command = "bowtie-build MT_genome/Danio_rerio.GRCz11.dna.chromosome.MT.fa " + index +  " 2> bowtie_build.txt"
   os.system(command)
   os.chdir("../")
   print("Downloading gtf annotation file .......")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/gtf/danio_rerio/Danio_rerio.GRCz11.92.chr.gtf.gz"
   print(command)
   os.system(command)
   os.system("gunzip -f Danio_rerio.GRCz11.92.chr.gtf.gz")
   infile= open("Danio_rerio.GRCz11.92.chr.gtf", "r")
   annotation = extract_MT(infile)


elif (args.species_name == "hsa") and (args.non_model==None):
   if not os.path.exists("human"):
    os.makedirs("human")
   species = "human"
   os.chdir("human")
   if not os.path.exists("MT_genome"):
    os.makedirs("MT_genome")
   os.chdir("MT_genome")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"
   print("Downloading Human Mitochondrial genome....................")
   print(command)
   os.system(command)
   os.system("gunzip -f *.gz")
   os.chdir("../")
   print("Creating bowtie index for Human mitochondrial genome")
   if not os.path.exists("bowtie-index"):
    os.makedirs("bowtie-index")
   index = "bowtie-index/human_MT_index"
   command = "bowtie-build MT_genome/Homo_sapiens.GRCh38.dna.chromosome.MT.fa " + index +  " 2> bowtie_build.txt"
   os.system(command)
   os.chdir("../")
   print("Downloading gtf annotation file .......")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.chr.gtf.gz"
   os.system(command)
   os.system("gunzip -f Homo_sapiens.GRCh38.92.chr.gtf.gz")
   infile= open("Homo_sapiens.GRCh38.92.chr.gtf", "r")
   annotation = extract_MT(infile)

elif (args.species_name == "mmu") and (args.non_model==None):
   if not os.path.exists("mouse"):
    os.makedirs("mouse")
   species = "mouse"
   os.chdir("mouse")
   if not os.path.exists("MT_genome"):
    os.makedirs("MT_genome")
   os.chdir("MT_genome")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.MT.fa.gz"
   print("Downloading Mouse Mitochondrial genome....................")
   print(command)
   os.system(command)
   os.system("gunzip -f *.gz")
   os.chdir("../")
   print("Creating bowtie index for Mouse mitochondrial genome")
   if not os.path.exists("bowtie-index"):
    os.makedirs("bowtie-index")
   index = "bowtie-index/mouse_MT_index"
   command = "bowtie-build MT_genome/Mus_musculus.GRCm38.dna.chromosome.MT.fa " + index +  " 2> bowtie_build.txt"
   os.system(command)
   os.chdir("../")
   print("Downloading gtf annotation file .......")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.chr.gtf.gz"
   os.system(command)
   os.system("gunzip -f Mus_musculus.GRCm38.92.chr.gtf.gz")
   infile= open("Mus_musculus.GRCm38.92.chr.gtf", "r")
   annotation = extract_MT(infile)


elif (args.species_name == "rno") and (args.non_model==None):
   if not os.path.exists("Rat"):
    os.makedirs("Rat")
   species = "rat"
   os.chdir("Rat")
   if not os.path.exists("MT_genome"):
    os.makedirs("MT_genome")
   os.chdir("MT_genome")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.chromosome.MT.fa.gz"
   print("Downloading Rat Mitochondrial genome....................")
   print(command)
   os.system(command)
   os.system("gunzip -f *.gz")
   os.chdir("../")
   print("Creating bowtie index for Rat mitochondrial genome")
   if not os.path.exists("bowtie-index"):
    os.makedirs("bowtie-index")
   index = "bowtie-index/Rat_MT_index"
   command = "bowtie-build MT_genome/Rattus_norvegicus.Rnor_6.0.dna.chromosome.MT.fa " + index +  " 2> bowtie_build.txt"
   os.system(command)
   os.chdir("../")
   print("Downloading gtf annotation file .......")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.92.chr.gtf.gz"
   os.system(command)
   os.system("gunzip -f Rattus_norvegicus.Rnor_6.0.92.chr.gtf.gz")
   infile= open("Rattus_norvegicus.Rnor_6.0.92.chr.gtf", "r")
   annotation = extract_MT(infile)

elif (args.species_name == "gga") and (args.non_model==None):
   if not os.path.exists("chicken"):
    os.makedirs("chicken")
   species = "chicken"
   os.chdir("chicken")
   if not os.path.exists("MT_genome"):
    os.makedirs("MT_genome")
   os.chdir("MT_genome")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/fasta/gallus_gallus/dna/Gallus_gallus.Gallus_gallus-5.0.dna.chromosome.MT.fa.gz"
   print("Downloading Rat Mitochondrial genome....................")
   print(command)
   os.system(command)
   os.system("gunzip -f *.gz")
   os.chdir("../")
   print("Creating bowtie index for Rat mitochondrial genome")
   if not os.path.exists("bowtie-index"):
    os.makedirs("bowtie-index")
   index = "bowtie-index/chicken_MT_index"
   command = "bowtie-build MT_genome/Gallus_gallus.Gallus_gallus-5.0.dna.chromosome.MT.fa " + index +  " 2> bowtie_build.txt"
   os.system(command)
   os.chdir("../")
   print("Downloading gtf annotation file .......")
   command = "wget ftp://ftp.ensembl.org/pub/release-92/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-5.0.92.chr.gtf.gz"
   os.system(command)
   os.system("gunzip -f Gallus_gallus.Gallus_gallus-5.0.92.chr.gtf.gz")
   infile= open("Gallus_gallus.Gallus_gallus-5.0.92.chr.gtf", "r")
   annotation = extract_MT(infile)

elif (args.species_name == "xen") and (args.non_model==None):
   if not os.path.exists("frog"):
    os.makedirs("frog")
   species = "frog"
   os.chdir("frog")
   if not os.path.exists("MT_genome"):
    os.makedirs("MT_genome")
   os.chdir("MT_genome")
   command = "wget ftp://ftp.xenbase.org/pub/Genomics/Sequences/Mitochondrial/Xlaevis_mito_seq.fa"
   print("Downloading Xenophus laevis Mitochondrial genome....................")
   print(command)
   os.system(command)
   os.chdir("../")
   print("Creating bowtie index for Xenophus laevis mitochondrial genome")
   if not os.path.exists("bowtie-index"):
    os.makedirs("bowtie-index")
   index = "bowtie-index/frog_MT_index"
   command = "bowtie-build MT_genome/Xlaevis_mito_seq.fa " + index +  " 2> bowtie_build.txt"
   os.system(command)
   os.chdir("../")
   print("Downloading gtf annotation file .......")
   command = "wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_GCF.gff3"
   os.system(command)
   infile= open("XENLA_9.2_GCF.gff3", "r")
   annotation = extract_MT(infile)

elif (args.species_name == "non_model"):
   species = "non_model"
   if not os.path.exists("non_model"):
    os.makedirs("non_model")
   os.chdir("non_model")
   if not os.path.exists("bowtie-index"):
    os.makedirs("bowtie-index")
   os.chdir("..")
   index = "bowtie-index/non_model_index"
   print("Creating bowtie index from the specified FASTA file")
   command = "bowtie-build " + args.FASTA + " " + "non_model/" + index +  " 2> bowtie_build.txt"
   os.system(command)
   infile= open(args.GTF, "r")
   annotation = extract_MT(infile)

else:
   sys.exit('\nERROR: Argument is not valid\n%s' %(docstring))

#extract filenames

if (args.files!=None) and (args.input_path!=None):
   sys.exit('\nERROR: --files and --input_path cannot be specified together\n%s'%(docstring))
elif (args.files!=None):
   files = arg.files
   files = [f for f in files if f.split("/")[-1].endswith(".fastq") or f.split("/")[-1].endswith(".fastq.gz")]
else:
   if (args.input_path==None):
      cwd = os.getcwd()
   else:
      cwd= args.input_path
   files = os.listdir(cwd)

if len(files) == 0:
   sys.exit('\nERROR: No FASTQ files in home directory\nEnsure that mtR_find is executed from the same place as the FASTQ files are present or specify the path to FASTQ files\n%s'%(docstring))

# s_RNA count step
if __name__ == "__main__":
 if args.RNA == "sRNA":
  if args.suspend == "no":
   pool = multiprocessing.Pool(processes = 30)
   result_list = pool.map(worker, [f for f in files if f.endswith(".fastq")])
   pool.close()
   pool.join()
  else:
   result_list = []
   files = [f for f in files if f.endswith(".fastq")]
   for f in files:
     result_list.append(worker(f))
 elif args.RNA == "lncRNA":
   result_list =[]
   file_list = [f for f in files if f.endswith(".fastq")]
   for f in file_list:
    infile= open(f)
    filename = str(f.split(".")[0].split("_")[0])
    fastq_lst = infile.readlines()[1::4]
    d = len(fastq_lst)
    fastq_lst = [line.strip() for line in fastq_lst]
    b = Counter(fastq_lst)
    result_list.append((d, filename, b))
 else:
   sys.exit('\nERROR: Argument is not valid\n%s' %(docstring))
 summed_counter = Counter()
 read_stat = []
 for x,y,z in result_list:
    print("Processing filename " + str(y))
    read_stat.append((y,x))
    summed_counter.update(z)
 read_sta = pd.DataFrame(list(read_stat), columns = ["filename", "total-count"])
 read_sta.to_csv("total_count_per_file.csv", sep = ",", index=False)
 summed_counter = dict(takewhile(lambda x: x[1] >= int(args.cutoff), summed_counter.most_common()))
 print("total number of unique sRNA sequences =  " + str(len(summed_counter)))


 summed_counter = pd.DataFrame.from_dict(summed_counter, orient='index').reset_index()
 summed_counter = summed_counter.rename(columns={'index':'read', 'total_count':'count'})

 summed_counter = pd.DataFrame(summed_counter)
 summed_counter.columns=['read', 'total_count']

 for x,y,z in result_list:
    df = pd.DataFrame.from_dict(z, orient='index').reset_index()
    df = df.rename(columns={'index':'read', 0:y})
    summed_counter = pd.merge(summed_counter, df, how='left', on='read')

 summed_counter = summed_counter.fillna(0)
 summed_counter.sort_values('read', inplace=True)
 summed_counter.to_csv("sRNA_count.csv", sep = ",", index=False)
 fasta_lst = summed_counter["read"].tolist()
 col_list=list(summed_counter)
 col_list[0]='sequence'
 summed_counter.columns=col_list

 ofile = open("master.fa", "w")
 for x in fasta_lst:
    ofile.write(">" +  str(x) +"\n" + str(x) + "\n")
 ofile.close()

 #map the sequences to mitochondrial genome index using bowtie
 command = "bowtie --best -v 1 -p 20 " + species + "/" + index + " -f master.fa -S filtered200_mt.sam 2> bowtie_log.txt"
 os.system(command)

 #seperate mapped and unmapped reads
 command = "samtools view -Sh -F 4 filtered200_mt.sam > filtered200_mt_mapped.sam"
 os.system(command)
 command = "samtools view -Sh -f 4 filtered200_mt.sam > filtered200_mt_unmapped.sam"
 os.system(command)

 #parsing SAM files
 command = "egrep -v '@HD|@SQ|@PG' filtered200_mt_mapped.sam > filtered200_mt_mapped.tsv"
 os.system(command)
 command = "egrep -v '@HD|@SQ|@PG' filtered200_mt_unmapped.sam > filtered200_mt_unmapped.tsv"
 os.system(command)




 lines = open("filtered200_mt_mapped.tsv", "r").readlines()
 allmtseq = pd.DataFrame([], columns=["sequence","type","bio-type","strand","orientation","annotation","Sequence alignment","Sequence start position(bp)","Sequence end position(bp)","gene-boundary:start(bp)","gene-boundary:end(bp)", "substitutions"])

 if (len(lines) != 0) and (args.RNA == "sRNA"):
   if args.suspend == "no":
    pool = multiprocessing.Pool()
    result_list = pool.map(partial(mt_annotator, kind = "normal_1_mismatch"), [line for line in lines])
    pool.close()
    pool.join()
   else:
    result_list = []
    for line in lines:
      result_list.append(mt_annotator(line,"normal_1_mismatch"))
   if len(result_list)==0:
      sys.exit('\nERROR: Problem with FASTQ or gtf file \n')
   allmtseq = pd.DataFrame(result_list, columns=["sequence","type","bio-type","strand","orientation","annotation","Sequence alignment","Sequence start position(bp)","Sequence end position(bp)","gene-boundary:start(bp)","gene-boundary:end(bp)", "substitutions"])
   print("Number of mtsRNAs mapping to mitochondrial genome with one-mismatch = " + str(len(allmtseq)))
   #mask CCA from unmapped reads
   lines = open("filtered200_mt_unmapped.tsv", "r").readlines()
   lines = [line.strip().split("\t")[0] for line in lines]
   selection=list(range(3,-4,-1))
   fo = open("CCA_fasta.fa","w")

   for k in lines:
    if (k[-3:] == 'CCA'):
       fo.write(">" + k[:-3] + "\n" + k[:-3] + "\n")
    else:
       continue
   if len(lines) > 0:
     #map the sequences to mitochondrial genome index using bowtie
     command = "bowtie --best -v 0 -p 20 " + species + "/" + index + " -f CCA_fasta.fa -S filtered200_CCA_mt.sam 2>> bowtie_log.txt"
     os.system(command)

     #seperate mapped and unmapped reads
     command = "samtools view -Sh -F 4 filtered200_CCA_mt.sam > filtered200_CCA_mt_mapped.sam"
     os.system(command)

     #parsing SAM files
     command = "egrep -v '@HD|@SQ|@PG' filtered200_CCA_mt_mapped.sam > filtered200_CCA_mt_mapped.tsv"
     os.system(command)

     lines2 = open("filtered200_CCA_mt_mapped.tsv", "r").readlines()
     if len(lines2) != 0:
        if args.suspend == "no":
         pool = multiprocessing.Pool()
         result_lst = pool.map(partial(mt_annotator, kind = "CCA_0_mismatch"), [line for line in lines2])
         pool.close()
         pool.join()
        else:
         result_lst = []
         for line in lines2:
           result_lst.append(mt_annotator(line,"CCA_0_mismatch"))
        selection=list(range(3,-4,-1))
        result_list = []
        for line in result_lst:
         if line[5] != "non-coding":
           if ((int(line[9])-int(line[7])) in selection) or ((int(line[10])-int(line[8])) in selection):
              line[0] = line[0] + "CCA"
              result_list.append(line)
           else:
              continue
        CCA = pd.DataFrame(result_list,columns=["sequence","type","bio-type", "strand","orientation","annotation","Sequence alignment","Sequence start position(bp)","Sequence end position(bp)","gene-boundary:start(bp)","gene-boundary:end(bp)", "substitutions"])
        if len(result_list) != 0:
           print("total number of CCA mitochondrial tRfs = " + str(len(CCA)))
           CCA["substitutions"] = "CCA"
           CCA.insert(1, "subtype", "")
           selection=list(range(3,-4,-1))
           for i in range(0,len(CCA),1):
            if len(CCA.at[i,'sequence']) - ((int(CCA.at[i,'gene-boundary:end(bp)'])-int(CCA.at[i,'gene-boundary:start(bp)']))/2) in selection:
             if (int(CCA.at[i,'gene-boundary:start(bp)'])-int(CCA.at[i,'Sequence start position(bp)'])) in selection:
                CCA.at[i,'subtype'] = 'tRNA-half-3'
             elif (int(CCA.at[i,'gene-boundary:end(bp)'])-int(CCA.at[i,'Sequence end position(bp)'])) in selection:
                CCA.at[i,'subtype']= 'tRNA-half-5'
            else:
             if (int(CCA.at[i,'gene-boundary:start(bp)'])-int(CCA.at[i,'Sequence start position(bp)'])) in selection:
                CCA.at[i,'subtype']= 'tRF-3'
             elif (int(CCA.at[i,'gene-boundary:end(bp)'])-int(CCA.at[i,'Sequence end position(bp)'])) in selection:
                CCA.at[i,'subtype']= 'tRF-5'
     else:
      CCA = pd.DataFrame([],columns=["sequence","type","bio-type", "strand","orientation","annotation","Sequence alignment","Sequence start position(bp)","Sequence end position(bp)","ge\
ne-boundary:start(bp)","gene-boundary:end(bp)", "substitutions"])
   else:
    CCA = pd.DataFrame([],columns=["sequence","type","bio-type", "strand","orientation","annotation","Sequence alignment","Sequence start position(bp)","Sequence end position(bp)","ge\
ne-boundary:start(bp)","gene-boundary:end(bp)", "substitutions"])

 if len(allmtseq) != 0:
  allmtseq.insert(1, "subtype", "")
  tRNA = allmtseq.loc[allmtseq["bio-type"]=="Mt_tRNA"]
  tRNA.reset_index(drop=True,inplace=True)


  for i in range(0,len(tRNA),1):
   if (tRNA.at[i,'strand'] == "H"):
    if len(tRNA.at[i,'sequence']) - ((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) in range(0,2):
     if (int(tRNA.at[i,'gene-boundary:start(bp)'])-int(tRNA.at[i,'Sequence start position(bp)'])) in selection:
        tRNA.at[i,'subtype']='tRNA-half-5'
     elif (int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'Sequence end position(bp)'])) in selection:
        tRNA.at[i,'subtype']='tRNA-half-3'
     elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) > int(tRNA.at[i,'Sequence start position(bp)']):
        tRNA.at[i,'subtype']='i-tRF-5'
     elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) <= int(tRNA.at[i,'Sequence start position(bp)']):
        tRNA.at[i,'subtype']='i-tRF-3'
    elif (int(tRNA.at[i,'gene-boundary:start(bp)'])-int(tRNA.at[i,'Sequence start position(bp)'])) in selection:
     tRNA.at[i,'subtype']='tRF-5'
    elif (int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'Sequence end position(bp)'])) in selection:
     tRNA.at[i,'subtype']='tRF-3'
    elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) > int(tRNA.at[i,'Sequence start position(bp)']):
     tRNA.at[i,'subtype']='i-tRF-5'
    elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) <= int(tRNA.at[i,'Sequence start position(bp)']):
     tRNA.at[i,'subtype']='i-tRF-3'
   elif tRNA.at[i,'strand'] == "L":
    if len(tRNA.at[i,'sequence']) - ((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) in range(0,2):
     if (int(tRNA.at[i,'gene-boundary:start(bp)'])-int(tRNA.at[i,'Sequence start position(bp)'])) in selection:
        tRNA.at[i,'subtype']='tRNA-half-3'
     elif (int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'Sequence end position(bp)'])) in selection:
        tRNA.at[i,'subtype']='tRNA-half-5'
     elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) > int(tRNA.at[i,'Sequence start position(bp)']):
        tRNA.at[i,'subtype']='i-tRF-5'
     elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) <= int(tRNA.at[i,'Sequence start position(bp)']):
        tRNA.at[i,'subtype']='i-tRF-3'
    elif (int(tRNA.at[i,'gene-boundary:start(bp)'])-int(tRNA.at[i,'Sequence start position(bp)'])) in selection:
     tRNA.at[i,'subtype']='tRF-3'
    elif (int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'Sequence end position(bp)'])) in selection:
     tRNA.at[i,'subtype']='tRF-5'
    elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) > int(tRNA.at[i,'Sequence start position(bp)']):
     tRNA.at[i,'subtype']='i-tRF-3'
    elif (((int(tRNA.at[i,'gene-boundary:end(bp)'])-int(tRNA.at[i,'gene-boundary:start(bp)']))/2) + int(tRNA.at[i,'gene-boundary:start(bp)'])) <= int(tRNA.at[i,'Sequence start position(bp)']):
     tRNA.at[i,'subtype']='i-tRF-5'


  na_count = allmtseq.loc[allmtseq['annotation'] == "non-coding"].copy()
  na_count["subtype"] = "non-coding"
  nontRNA = allmtseq.loc[(allmtseq["bio-type"] == "protein_coding")|(allmtseq["bio-type"] == "Mt_rRNA")]
  nontRNA.reset_index(drop=True,inplace=True)

  for i in range(0,len(nontRNA),1):
   if (nontRNA.at[i,'strand'] == "H"):
    if len(nontRNA.at[i,'sequence']) - ((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) in range(0,2):
     if (int(nontRNA.at[i,'gene-boundary:start(bp)'])-int(nontRNA.at[i,'Sequence start position(bp)'])) in selection:
        nontRNA.at[i,'subtype']="5'-half"
     elif (int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'Sequence end position(bp)'])) in selection:
        nontRNA.at[i,'subtype']="3'-half"
     elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) > int(nontRNA.at[i,'Sequence start position(bp)']):
        nontRNA.at[i,'subtype']='i-5prime'
     elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) <= int(nontRNA.at[i,'Sequence start position(bp)']):
        nontRNA.at[i,'subtype']='i-3prime'
    elif (int(nontRNA.at[i,'gene-boundary:start(bp)'])-int(nontRNA.at[i,'Sequence start position(bp)'])) in selection:
     nontRNA.at[i,'subtype']="5prime"
    elif (int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'Sequence end position(bp)'])) in selection:
     nontRNA.at[i,'subtype']="3prime"
    elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) > int(nontRNA.at[i,'Sequence start position(bp)']):
     nontRNA.at[i,'subtype']='i-5prime'
    elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) <= int(nontRNA.at[i,'Sequence start position(bp)']):
     nontRNA.at[i,'subtype']='i-3prime'
   elif (nontRNA.at[i,'strand'] == "L"):
    if len(nontRNA.at[i,'sequence']) - ((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) in range(0,2):
     if (int(nontRNA.at[i,'gene-boundary:start(bp)'])-int(nontRNA.at[i,'Sequence start position(bp)'])) in selection:
        nontRNA.at[i,'subtype']="3'-half"
     elif (int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'Sequence end position(bp)'])) in selection:
        nontRNA.at[i,'subtype']="5'-half"
     elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) > int(nontRNA.at[i,'Sequence start position(bp)']):
        nontRNA.at[i,'subtype']='i-5prime'
     elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) <= int(nontRNA.at[i,'Sequence start position(bp)']):
        nontRNA.at[i,'subtype']='i-3prime'
    elif (int(nontRNA.at[i,'gene-boundary:start(bp)'])-int(nontRNA.at[i,'Sequence start position(bp)'])) in selection:
     nontRNA.at[i,'subtype']="3prime"
    elif (int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'Sequence end position(bp)'])) in selection:
     nontRNA.at[i,'subtype']="5prime"
    elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) > int(nontRNA.at[i,'Sequence start position(bp)']):
     nontRNA.at[i,'subtype']='i-3prime'
    elif (((int(nontRNA.at[i,'gene-boundary:end(bp)'])-int(nontRNA.at[i,'gene-boundary:start(bp)']))/2) + int(nontRNA.at[i,'gene-boundary:start(bp)'])) <= int(nontRNA.at[i,'Sequence start position(bp)']):
     nontRNA.at[i,'subtype']='i-5prime'


  filter_lst = []
  for i in range(0,len(annotation)):
    if ((annotation[i][4] == "Mt_tRNA") and (annotation[i][2] == "+")) and (annotation[i][2] == annotation[i+1][2]):
     filter_lst.append((annotation[i][3],annotation[i+1][3]))

  for x,v in filter_lst:
   for i in range(0,len(tRNA),1):
    if (tRNA.at[i,'subtype'] == 'tRF-5') and (tRNA.at[i,'annotation'] == v) and (tRNA.at[i,'strand'] == "H"):
     z = "tRF-1"
     tRNA.at[i,'subtype']=z
     tRNA.at[i,'annotation'] = str(x)
    if (tRNA.at[i,'subtype'] == "5'-half") and (tRNA.at[i,'annotation'] == v) and (tRNA.at[i,'strand'] == "H"):
     z = "tRF-1"
     tRNA.at[i,'subtype']=z
     tRNA.at[i,'annotation']=str(x)

  for x,v in filter_lst:
   for i in range(0,len(nontRNA),1):
    if (nontRNA.at[i,'subtype'] == '5prime') and (nontRNA.at[i,'annotation'] == v) and (nontRNA.at[i,'strand'] == "H"):
     z = "tRF-1"
     nontRNA.at[i,'subtype']=z
     nontRNA.at[i,'annotation']=str(x)
  foy_lst = [CCA, nontRNA, na_count]
  for x in foy_lst:
   if len(y) == 0:
      foy_lst.remove(x)
  allmtseq = tRNA.append(foy_lst)[tRNA.columns.tolist()]
  allmtseq=pd.merge(allmtseq, summed_counter, how='left', on='sequence')
  allmtseq.insert(1, "Specific-ID", "")
  allmtseq.insert(2, "General-ID", "")
  allmtseq.insert(3, "temp_c", "")
  allmtseq.insert(4, "temp_d", "")
  allmtseq.insert(5, "temp_e", "")
  allmtseq.insert(6, "temp_f", "")
  fnames = {"anti-sense":"as","sense": ""}
  gnames = {"normal_1_mismatch":"", "CCA_0_mismatch": "-CCA"}
  nnames = {'tRNA-half-5': 'tRH-5', 'tRNA-half-3': 'tRH-3', 'i-tRF-5': 'i-tRF-5', 'i-tRF-3': 'i-tRF-3', 'tRF-3': 'tRF-3', 'tRF-5': 'tRF-5', 'tRF-1': 'tRF-1', 'i-3prime':'i-3p', 'i-5prime': 'i-5p', '5prime': '5p', '3prime': '3p', "5'-half": '5H', "3'-half": '3H', "non-coding": 'nc' }  
  allmtseq["temp_c"] = allmtseq["orientation"].map(fnames)
  allmtseq["temp_d"] = allmtseq["type"].map(gnames)
  allmtseq["temp_e"] = allmtseq.apply(lambda x: len(x["sequence"]), axis =1)
  allmtseq = allmtseq.loc[allmtseq["temp_e"]<51]
  allmtseq["temp_e"] = allmtseq.astype(str)
  allmtseq["temp_f"] = allmtseq["subtype"].map(nnames)
  allmtseq["Specific-ID"]= args.species_name + "|" + "mt-sRNA" + "|" + allmtseq["annotation"]+ "|" + allmtseq["strand"] + "|" + allmtseq["temp_c"] + allmtseq["Sequence start position(bp)"] + "|" +  allmtseq["temp_f"] + "|" + allmtseq["temp_e"] + "|" +allmtseq["substitutions"] 
  allmtseq["General-ID"]= args.species_name + "|" + "mt-sRNA" + "|" + allmtseq["annotation"]+ "|" + allmtseq["strand"] + "|" + allmtseq["temp_c"] + allmtseq["Sequence start position(bp)"] + allmtseq["temp_d"] 
  allmtseq.drop(["temp_c","temp_d","temp_e","temp_f"], inplace = True, axis =1)
  if len(allmtseq) > 0:
   allmtseq.to_csv("mt-sRNA_count.csv", sep = ",", index=False)
  print("Total number of mtsRNAs is " + str(len(allmtseq)))


 #long non coding RNA section

 if (len(lines) != 0) and (args.RNA == "lncRNA"):
   pool = multiprocessing.Pool()
   result_list = pool.map(partial(mt_annotator2, kind = "mt-lncRNA"), [line for line in lines])
   pool.close()
   pool.join()
   allmtseq2 = pd.DataFrame(result_list, columns=["sequence","type","bio-type","strand","orientation","annotation","Sequence alignment","Sequence start position(bp)","Sequence end position(bp)","gene-boundary:start(bp)","gene-boundary:end(bp)", "substitutions"])
   allmtseq2=pd.merge(allmtseq2, summed_counter, how='left', on='sequence')
   allmtseq2.insert(1, "Specific-ID", "")
   allmtseq2.insert(2, "Length", "")
   allmtseq2.insert(2, "General-ID", "")
   allmtseq2["Length"] = allmtseq2.apply(lambda x: len(x["sequence"]), axis = 1)
   allmtseq2["Length"] = allmtseq2["Length"].astype(int)
   allmtseq2 = allmtseq2.loc[allmtseq2["Length"]>50]
   allmtseq2 = allmtseq2.loc[allmtseq2["Length"]> int(args.filter)]
   allmtseq2["Length"] = allmtseq2["Length"].astype(str)
   allmtseq2["Specific-ID"]= args.species_name + "|" + allmtseq2["type"] + "|" + allmtseq2["annotation"]+ "|" + allmtseq2["strand"] + "|" + allmtseq2["Sequence start position(bp)"] + "|" + allmtseq2["Length"].apply(str) +  allmtseq2["substitutions"]
   allmtseq2["General-ID"]= args.species_name + "|" + allmtseq2["type"] + "|"  + allmtseq2["annotation"]+ "|" + allmtseq2["strand"] + "|" + allmtseq2["Sequence start position(bp)"] 
   allmtseq2.sort_values('sequence', inplace=True)
   allmtseq2.reset_index(drop=True,inplace=True)
   allmtseq2.to_csv("lnc_RNA.csv", sep = ",", index=False)
   sel_list = ["sequence","total_count"]
   sel_list.extend([f.split(".")[0] for f in files if f.endswith(".fastq")])
   csv = allmtseq2.loc[:,sel_list]
   asan = pd.DataFrame([],columns=sel_list)

   #initialize i=0
   i=0

   while(i<=(len(csv)-2)):
      seed = csv.at[i,'sequence']
      match = csv.at[i+1,'sequence']
      k = len(match)-len(seed)
      if seed != match[:-k]:
         nasa=csv.iloc[i,]
         asan=asan.append(nasa, ignore_index=True)
         i=i+1
      if seed == match[:-k]:
         j=i+1
         while(j<=(len(csv))-1):
          seed = csv.at[i,'sequence']
          match = csv.at[j,'sequence']
          m = len(match)-len(seed)
          if seed == match[:-m]:
             j = j+1
          else:
             break
         nasa=csv.iloc[i:j,1:32]
         tot = nasa.sum(axis=0)
         asan=asan.append(tot, ignore_index=True)
         asan['sequence'] = asan['sequence'].astype('str')
         seq = csv.at[i,'sequence']
         z= len(asan)-1
         asan.at[z,'sequence']=seq
         i = j
   allmtseq2=allmtseq2.iloc[:,0:13]
   allmtseq3=allmtseq2.loc[allmtseq2["sequence"].isin(asan["sequence"].tolist())]
   allmtseq3=pd.merge(allmtseq3, asan, how='left', on='sequence')
   allmtseq3.to_csv("lnc_RNA1.csv", sep = ",", index=False)
   print("Total number of mtlncRNA is " + str(len(allmtseq3)))

 if (args.graphical_output == "yes"):
  if args.RNA == "lncRNA":
     master = allmtseq3
  else:
     master = allmtseq
  metadata = pd.read_csv(args.metadata, sep = None)
  condition = list(metadata)[1]
  df=metadata.groupby(condition)

  d = {}
  colors = ['red', 'green', 'navy','yellow','orange','black','brown','salmon','steelblue','maroon','cyan','deeppink','darkslategray','wheat','rosybrown','olive']
  if len(df) > len(colors):
   sys.exit('\nERROR: NUmber of condition excees the number on the pre-defined color palette. Please enter colors manually using --colors option\n')
  for name,ind_df in df:
     d[name] = ind_df[list(ind_df)[0]].tolist()

  new_lst = ['strand', 'orientation', 'annotation']


  for k,v in list(d.items()):
     merged_list = []
     v = [ item.split(".")[0] for item in v]
     merged_list = new_lst + v
     subset_df = master.loc[:,merged_list]
     subset_df["total_count"]= subset_df[v].sum(axis=1)
     df_strand = subset_df.groupby("strand").sum()
     df_strand.reset_index(level=0, inplace=True)
     output = k + "strand.csv"
     df_strand.to_csv(output, sep = ",", index=False)
     df_orientation = subset_df.groupby("orientation").sum()
     df_orientation.reset_index(level=0, inplace=True)
     output = k + "orientation.csv"
     df_orientation.to_csv(output, sep = ",", index=False)
     df_annotation = subset_df.groupby("annotation").sum()
     df_annotation.reset_index(level=0, inplace=True)
     output = k + "annotation.csv"
     df_annotation.to_csv(output, sep = ",", index=False)
     venn_dict = dict(list(zip(df_strand.strand, df_strand.total_count)))
     plt.rcParams['font.size'] = 9
     labels2 = ['Heavy Strand(%)', 'Light Strand(%)']
     sizes4 = [venn_dict["H"], venn_dict["L"]]
     fig1, ax1 = plt.subplots()
     explode1 = (0, 0.1)
     ax1.pie(sizes4, explode=explode1, autopct='%1.1f',shadow=False, startangle=90, colors=['limegreen','red'])
     ax1.axis('equal')
     ax1.set_title(k)
     ax1.legend(labels2,fontsize=8, loc='best')
     plt.tight_layout()
     plt.savefig(str(k) + '_strand.png', dpi = 300, bbox_inches='tight')
     plt.close()
     opacity = 0.4
     color = random.choice(colors)
     colors.remove(color)
     rects1 = df_annotation.plot.bar(x="annotation", y="total_count", rot = 90,legend = False,color=color, fontsize = 7)
     rects1.set_xlabel('MtDNA genes',fontsize=12)
     rects1.set_ylabel('Number of reads',fontsize=12)
     rects1.set_title("Distribution of reads", fontsize = 14)
     plt.tight_layout()
     praph = str(output.split(".")[0]) + ".png"
     plt.savefig(praph, dpi = 300, bbox_inches='tight')
     plt.close()
