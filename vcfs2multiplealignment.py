#Import Libraries

import sys
import os
import argparse
import gzip



parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-ref', action="store", type=str, required=True, dest='referencefile', default="", help='Reference file')
parser.add_argument('-vcfs', action="store", type=str, required=True, dest='vcfs', nargs="+", default="", help='vcf')
parser.add_argument('-maj_s', action="store", type=float, dest='maj_s', default=0.8, help='Support for accepting majority variant')
parser.add_argument('-d', action="store", type=int, dest='depth', default=50, help='Depth threshold for including a position')
parser.add_argument('-min_s', action="store", type=float, dest='min_s', default=0.20, help='Support for accepting minority variant')
parser.add_argument('-gap_s', action="store", type=float, dest='gap_s', default=0.70, help='Support for accepting gaps')

args = parser.parse_args()

referencefile = args.referencefile
vcfs = args.vcfs
depth = args.depth
maj_s = args.maj_s
min_s = args.min_s
gap_s = args.gap_s

 # SHOULD BE MADE TO HANDLE MULTIPLE BASES PER VCF POSITION!!

def main():
    sequencelist = []
    headerlist = []

    for file in vcfs:
        header, refsequence = loadSequence(referencefile)  # Return header as string, sequence as list
        header, sequence = makeSingleSequnece(header, refsequence, file)
        headerlist.append(header)
        sequencelist.append(sequence)


    positions = len(sequence)
    sequencelist_len = len(sequencelist)
    for i in range(positions):
        max = 0
        for t in range(sequencelist_len):
            if len(sequencelist[t][i]) > max:
                max = len(sequencelist[t][i])
        if max > 1: #insertion has occured
            for t in range(sequencelist_len):
                if max > len(sequencelist[t][i]): #Select positions needing gaps
                    diff = max - len(sequencelist[t][i])
                    gaps = "-" * diff
                    sequencelist[t][i] = sequencelist[t][i] + gaps
    for i in range(len(sequencelist)):
        consensus = "".join(sequencelist[i])
        print (headerlist[i])
        print (consensus)







def makeSingleSequnece(header, sequence, vcf):
    vcfList , vcfName = loadVCF(vcf)
    header, sequence = consensusMaker(header, sequence, vcfList, depth, maj_s, min_s, vcfName, gap_s)
    return header, sequence


def consensusMaker(header, sequence, vcfList, depth, maj_s, min_s, vcfName, gap_s):
    for position in vcfList:
        #print (position)
        if position[3] == "<->":
            position[3] = "-"
        if position[4] == "<->":
            position[4] = "-"
        vcfInfo = position[7].split(";")
        if float(vcfInfo[0][3:]) >= depth: #Check depth
            #Check for minority variant:
            variants = vcfInfo[5][4:]
            vcfInfo = position[7].split(";")
            dp = int(vcfInfo[0][3:])
            variantCount = variants.split(",")
            gaps_support = int(variantCount[-1]) / dp
            sort_variantCount = []
            for i in range(len(variantCount)):
                sort_variantCount.append(int(variantCount[i]))
            sort_variantCount.sort()
            minority_index = variantCount.index(str(sort_variantCount[-2]))
            dp = int(vcfInfo[0][3:])
            variantList = ['A', 'C', 'G', 'T', 'N', '-']
            minorityVariant = variantList[minority_index]
            minority_depth = int(sort_variantCount[-2]) / dp
            if position[3] != "-" or position[4] != "-": #spørg PLAN om det her
                if minority_depth > min_s: #Minority Variant found
                    if position[1] != "0":
                        if position[4] == "-":  # Handle deletion
                            sequence[int(position[1])-1] = "-"
                        elif gaps_support >= gap_s:
                            sequence[int(position[1])-1] = "-"
                        else: #minority variant
                            sequence[int(position[1])-1] = minorityVariant
                    else: #insertion
                        sequence[int(previousPositions[1])-1] = sequence[int(previousPositions[1])-1] + minorityVariant
                elif float(vcfInfo[2][3:]) >= maj_s: #Acceptable support for majority variants
                    if position[4] == "-":
                        position[4] = "-"
                    if position[1] != "0":
                        if position[4] == "-":  # Handle deletion
                            sequence[int(position[1])-1] = "-"
                        elif gaps_support >= gap_s:
                            sequence[int(position[1])-1] = "-"
                        else:  # Majority call
                            sequence[int(position[1])-1] = position[4]
                    else:  # insertion
                        if position[3] != "-" or position[4] != "-":
                            sequence[int(previousPositions[1])-1] = sequence[int(previousPositions[1])-1] + position[4]
        if position[1] != "0":
            previousPositions = position
    header = ">" + vcfName
    return header, sequence

def loadSequence(referencefile):
    sequence = []
    if referencefile[-2:] == 'gz':
        infile = gzip.open(referencefile, 'rb')
        for line in infile:
            line = line.decode()
            line = line.rstrip()
            if line[0] == ">":
                header = line
            else:
                for nucleotide in line:
                    sequence.append(nucleotide)
        infile.close()

    else:
        infile = open(referencefile, 'r')
        for line in infile:
            line = line.rstrip()
            if line[0] == ">":
                header = line
            else:
                for nucleotide in line:
                    sequence.append(nucleotide)
        infile.close()
    return header, sequence

def loadVCF(vcf):
    vcfList = []
    if vcf[-2:] == 'gz':
        infile = gzip.open(vcf, 'rb')
        for line in infile:
            line = line.decode()
            line = line.rstrip()
            if line[0] != "#":
                line = line.split("\t")
                vcfList.append(line)
    else:
        infile = open(vcf, 'r')
        for line in infile:
            line = line.rstrip()
            if line[0] != "#":
                line = line.split("\t")
                vcfList.append(line)
    if "/" in vcf:
        vcfName = vcf.split("/")[-1]
    else:
        vcfName = vcf
    return vcfList, vcfName

def checkInput(variantCode):
    if variantCode > 3 or variantCode < 1:
        sys.exit("VariantCode must be either 1, 2 or 3")



if __name__ == '__main__':
    main()