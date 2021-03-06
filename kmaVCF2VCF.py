

import sys
import os
import argparse
import gzip

#Formats KMA's VCF from AD6 to DP4 and fixes KMA indel format.

parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-vcf', action="store", type=str, required=True, dest='vcf', default="", help='kma vcf')
parser.add_argument('-d', action="store", type=int, dest='depth', default=100, help='Depth threshold for including a position')
#parser.add_argument('-maj_s', action="store", type=float, dest='maj_s', default=0.7, help='Support for accepting majority variant')
parser.add_argument('-min_s', action="store", type=float, dest='min_s', default=1, help='Support for accepting minority variant. To call minority variants, set min_s to ~0.2')
parser.add_argument('-gap_s', action="store", type=float, dest='gap_s', default=0.40, help='Support for accepting gaps')

args = parser.parse_args()

vcfheader = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n##INFO=<ID=AD6,Number=6,Type=Integer,Description=\"Count of all alternative alleles: A,C,G,T,N,-\">\">\n##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description=\"Indicates that the variant is a consensus variant (as opposed to a low frequency variant).\">\n##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\"Homopolymer length to the right of report indel position\">"
vcf = args.vcf
depth = args.depth
min_s = args.min_s
gap_s = args.gap_s


def main():
    vcflist = loadVCF(vcf)
    vcflist = convertVCF(vcflist, min_s)
    vcflist = indelTag(vcflist)
    print (vcfheader)
    for i in range(len(vcflist)):
        vcflist[i][3] = vcflist[i][3].upper()
        vcflist[i][4] = vcflist[i][4].upper()
        if vcflist[i][3] != vcflist[i][4]:
            print ("\t".join(vcflist[i]))

def indelTag(vcflist):
    newVCFlist = []
    for position in vcflist:
        if len(position[3]) > 1 or len(position[4]) > 1:
            position[7] = position[7] + ";INDEL"
        newVCFlist.append(position)
    return newVCFlist

def convertVCF(vcflist, min_s):
    newVCFlist = []
    for position in vcflist:
        vcfInfo = position[7].split(";")
        if float(vcfInfo[0][3:]) >= depth:
            positionType, positionvcfInfo = indentifyPositionType(position, gap_s)
            minorityVariant, minority_depth = calculateMinor(vcfInfo)
            newVCFlist = handlePosition(position, positionType, minorityVariant, minority_depth, min_s, newVCFlist, positionvcfInfo)
    return newVCFlist

def handlePosition(position, positionType, minorityVariant, minority_depth, min_s, newVCFlist, positionvcfInfo):
    indelFlag = False

    if positionType == "variant_majority":
        if minority_depth >= min_s:
            if minorityVariant != "-":
                position[4] = minorityVariant
                newVCFlist.append(position)
            else:
                if newVCFlist == []:
                    newVCFlist.append(position)
                else:
                    newVCFlist[-1][4] += position[4]
        else:
            newVCFlist.append(position)
    elif positionType == "deletion_majority":
        newVCFlist[-1][3] += position[3]
        newVCFlist[-1][7] = positionvcfInfo
    elif positionType == "insertion_majority":
        newVCFlist[-1][4] += position[4]
        newVCFlist[-1][7] = positionvcfInfo
    elif positionType == "minor_insertion":
        pass
    else:
        newVCFlist.append(position)

    return newVCFlist

def calculateMinor(vcfInfo):
    variants = vcfInfo[5][4:]
    variantCount = variants.split(",")
    sort_variantCount = []
    for i in range(len(variantCount)):
        sort_variantCount.append(int(variantCount[i]))
    sort_variantCount.sort()
    minority_index = variantCount.index(str(sort_variantCount[-2]))
    dp = int(vcfInfo[0][3:])
    variantList = ['A', 'C', 'G', 'T', 'N', '-']
    minorityVariant = variantList[minority_index]
    minority_depth = int(sort_variantCount[-2]) / dp
    return minorityVariant, minority_depth

def indentifyPositionType(position, gap_s):
    vcfInfo = position[7].split(";")
    dp = int(vcfInfo[0][3:])
    ad6 = vcfInfo[5][4:]
    ad6list = ad6.split(",")
    if position[1] != "0":
        gaps_support = int(ad6list[-1])/dp
        #print (gaps_support)
        if position[4] == "<->":  # deletion:
            positionType = "deletion_majority"
        elif gaps_support >= gap_s:
            positionType = "deletion_majority"
        else:
            positionType = "variant_majority"
    else:
        if position[3] == "<->" and position[4] != "<->":  # Normal Insertion majority
            positionType = "insertion_majority"

        else:
            positionType = "minor_insertion"
    vcfInfo = position[7]
    return positionType, vcfInfo

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

