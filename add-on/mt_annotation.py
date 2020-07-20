#!/usr/bin/env python
from builtins import range
import argparse

#extract mtDNA annotation from gtf file
parser = argparse.ArgumentParser(usage = "\n" +"python %(prog)s <input_file> <output_file> \n" +"\n" + "Usage examples\
:\npython mt_annotation.py mouse.gtf mouse_mt.gtf n\n"+ "\n" + "Description:\nExtracts mtDNa annotaions from GTF file\n", add_help=False)
required = parser.add_argument_group('required arguments')
required.add_argument("input_file", help = "path to input GTF file")
required.add_argument("output_file", help = "path to output GTF file")
optional = parser.add_argument_group('optional arguments')
optional.add_argument("-h", "--help", action='help', help='print help message')
args = parser.parse_args()

#read input file
infile = open(args.input_file, "r")
#skip the header lines
lines = infile.readlines()[5:]

# and output file name in the second
outfile = args.output_file
fo = open(outfile,"w")
fo.write("#!genome-build GRCz10"+ "\n"+ "#!genome-version GRCz10"+ "\n"+"#!genome-date 2014-09"+ "\n"+"#!genome-build-accession NCBI:GCA_000002035.3"+ "\n"+"#!genebuild-last-updated 2016-07"+"\n")
fo.write("seqname" + "\t"+ "source" + "\t" + "feature" +"\t" + "start" + "\t"+ "end" + "\t" + "score" +"\t" + "\t"+ "strand" + "\t" + "frame" +"\t"+ "attribute" + "\t" + "\n" )
for line in lines:
 line= line.strip().split("\t")
 if (line[0] == 'MT' and line[2]== 'gene'):
     #print line
     ten =''
     for i in range(0,len(line)):
       ten=ten+ line[i] +"\t"
     #print ten  
     fo.write(ten+"\n")
fo.close()
infile.close()
