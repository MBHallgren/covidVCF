#Import Libraries

import sys
import os
import argparse
import gzip
import json

#Identifies variants according to inputput type, and can give a binary output when searching for a specific variant

parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-i', action="store", type=str, required=True, dest='input', default="", help='input file')

args = parser.parse_args()

#Initialize
input = args.input
output_dir = "/Users/malhal/dev/data/outdir/" #Set
script_path = "/Users/malhal/dev/covidVCF/"
kma_db_path = "" #Set

def main(input):
    identifyReadType(input)

def identifyReadType(input):
    # Get inputfiles splitted in Nanopore and non-nanopre sequences.
    cmd = "%s{}fingerseq -i %s | grep -v \"^#\"> %s/fingerReport.tsv".format(script_path) % (
        paths['scripts'], " ".join(inputFiles), paths['outputs']);
    os.system(cmd);
    infile = open(paths['outputs'] + "/fingerReport.tsv", "r");
    Illumina_files = [];
    Nano_files = [];
    Asm_files = [];
    PE = 0;
    fsaInput = False;
    for line in infile:
        line = line.rstrip();
        info = line.split("\t");
        if (info[1] == "Nanopore"):
            Nano_files.append(info[0]);
        elif (info[3] != "Na"):
            Illumina_files.append(info[0]);
            PE += 1;
        elif (info[1] == "fastA"):
            Asm_files.append(info[0]);
            fsaInput = True;
        else:
            Illumina_files.append(info[0]);

    print (Illumina_files)
    print (Nano_files)

def kmaRun(input):
    cmd = "{} -i {} -Mt1 1 -bc 1.0 -bcd 100000 -1t1 -mrs 0.25 -vcf -nf -bcNano ".format()

def errorHandle(input):
    return True

if __name__ == '__main__':
    main(input)