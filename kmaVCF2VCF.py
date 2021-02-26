#Import Libraries

import sys
import os
import argparse
import gzip

#Formats KMA's VCF from AD6 to DP4 and fixes KMA indel format.

parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-vcf', action="store", type=str, required=True, dest='vcf', default="", help='kma vcf')


args = parser.parse_args()

vcfheader = "##reference=/home/share/pkrisz5-sarscov2veo/data/ref/sars2/NC_045512.2.fa\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description=\"Indicates that the variant is a consensus variant (as opposed to a low frequency variant).\">\n##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\"Homopolymer length to the right of report indel position\">"
basepairDict = dict()
basepairDict['A'] = 'T'
basepairDict['T'] = 'A'
basepairDict['C'] = 'G'
basepairDict['G'] = 'C'
vcf = args.vcf

def main():
    vcflist = loadVCF(vcf)
    vcflist = convertVCF(vcflist)
    #vcflist = covertBaseCalls(vcflist)
    print (vcfheader)
    for i in range(len(vcflist)):
        print ("\t".join(vcflist[i]))

def covertBaseCalls(vcflist):
    newVCFlist = []
    for position in vcflist:
        dp4 = [0,0,0,0]
        ad6 = position[7].split(";")[5][4:].split(",")
        if position[3] == "A":
            dp4[0] = int(ad6[0])/2
            dp4[1] = int(ad6[0])/2
        elif position[3] == "C":
            dp4[0] = int(ad6[1])/2
            dp4[1] = int(ad6[1])/2
        elif position[3] == "G":
            dp4[0] = int(ad6[2])/2
            dp4[1] = int(ad6[2])/2
        elif position[3] == "T":
            dp4[0] = int(ad6[3])/2
            dp4[1] = int(ad6[3])/2
        elif position[3] == "N":
            dp4[0] = int(ad6[4])/2
            dp4[1] = int(ad6[4])/2
        elif position[3] == "-":
            dp4[0] = int(ad6[5])/2
            dp4[1] = int(ad6[5])/2

        if position[4].upper() == "A":
            dp4[2] = int(ad6[0])/2
            dp4[3] = int(ad6[0])/2
        elif position[4].upper() == "C":
            dp4[2] = int(ad6[1])/2
            dp4[3] = int(ad6[1])/2
        elif position[4].upper() == "G":
            dp4[2] = int(ad6[2])/2
            dp4[3] = int(ad6[2])/2
        elif position[4].upper() == "T":
            dp4[2] = int(ad6[3])/2
            dp4[3] = int(ad6[3])/2
        elif position[4].upper() == "N":
            dp4[2] = int(ad6[4])/2
            dp4[3] = int(ad6[4])/2
        elif position[4].upper() == "-":
            dp4[2] = int(ad6[5])/2
            dp4[3] = int(ad6[5])/2
        for i in range(len(dp4)):
            dp4[i] = str(dp4[i])
        vcfinfo = position[7].split(";")
        vcfinfo[5] = "DP4=" + ",".join(dp4)
        position[7] = ";".join(vcfinfo)
        newVCFlist.append(position)
    return newVCFlist


def convertVCF(vcflist):
    newVCFlist = []
    for position in vcflist:
        if position[1] != "0":
            if position[4] == "<->": #deletion:
                if position[3] != "<->":
                    newVCFlist[-1][3] += position[3]
            else:
                newVCFlist.append(position)
        else:
            if position[3] == "<->" and position[4] != "<->": #Insertion
                newVCFlist[-1][4] += position[4]

    return newVCFlist



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
    return vcfList

if __name__ == '__main__':
    main()