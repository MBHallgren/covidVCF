

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
parser.add_argument('-gap_s', action="store", type=float, dest='gap_s', default=0.70, help='Support for accepting gaps')

args = parser.parse_args()

vcfheader = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n##INFO=<ID=DP2,Number=1,Type=Integer,Description=\"Count of the alternative allele\">\">\n##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description=\"Indicates that the variant is a consensus variant (as opposed to a low frequency variant).\">\n##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\"Homopolymer length to the right of report indel position\">"
vcf = args.vcf
depth = args.depth
min_s = args.min_s
gap_s = args.gap_s


def main():
    vcflist = loadVCF(vcf)
    vcflist = removeDoubleGaps(vcflist, depth)
    vcflist = insertMinor(vcflist, min_s)
    vcflist = convertVCF(vcflist, min_s)
    vcflist = indelTag(vcflist)
    print (vcfheader)
    for i in range(len(vcflist)):
        vcflist[i][3] = vcflist[i][3].upper()
        vcflist[i][4] = vcflist[i][4].upper()
        if vcflist[i][3] != vcflist[i][4]:
            print ("\t".join(vcflist[i]))

def insertMinor(vcflist, min_s):
    newvcflist = []
    for position in vcflist:
        newvcflist.append(position)
        vcfInfo = position[7].split(";")
        if vcfInfo[-1][0:3]!= "DP2":
            minorityVariant, minority_depth, minority_count = calculateMinor(vcfInfo)
            if minority_depth >= min_s:
                if minorityVariant == "-":
                    position[4] = "<->"
                    newvcflist.append(position)
                else:
                    if minorityVariant.upper() == position[3]:
                        pass
                    else:
                        position[4] = minorityVariant
                        dp = int(vcfInfo[0][3:])
                        vcfInfo[-1] = "DP2:" + str(minority_count)
                        vcfInfo[1] = "AD=" + str(minority_count)
                        vcfInfo[2] = "AF=" + str(minority_count / dp)[0:4]
                        vcfInfo = ";".join(vcfInfo)
                        position[7] = vcfInfo
                        position[5] = "0"
                        newvcflist.append(position)
    return vcflist

def removeDoubleGaps(vcflist, depth):
    newvcflist = []
    for position in vcflist:
        vcfInfo = position[7].split(";")
        dp = int(vcfInfo[0][3:])
        if position[3] == "<->" and position[4] == "<->":
            pass
        elif depth > dp: #Exclude non supported positions
            pass
        else:
            newvcflist.append(position)
    return newvcflist

def indelTag(vcflist):
    newVCFlist = []
    for position in vcflist:
        if len(position[3]) > 1 or len(position[4]) > 1:
            position[7] = position[7] + ";INDEL"
        newVCFlist.append(position)
    return newVCFlist

def checkIndel(position, gap_s, minority_depth, min_s, minorityVariant):
    indelCheck = False
    type = None
    if position[3] == "<->":
        type = "insertion"
        indelCheck = True
    if position[4] == "<->":
        type = "deletion"
        indelCheck = True
    vcfInfo = position[7].split(";")
    dp = int(vcfInfo[0][3:])
    ad6 = vcfInfo[5][4:]
    ad6list = ad6.split(",")
    gaps_support = int(ad6list[-1]) / dp
    if gaps_support > gap_s:
        type = "deletion"
        indelCheck = True
    if minority_depth >= min_s:
        if minorityVariant == "-":
            type = "deletion"
            indelCheck = True
    return indelCheck, type

def convertVCF(vcflist, min_s):
    newVCFlist = []
    indellist = []
    prev_position = None
    for position in vcflist:
        vcfInfo = position[7].split(";")
        if vcfInfo[-1][0:3]== "DP2":
            indelCheck = False
        else:
            minorityVariant, minority_depth, minority_count = calculateMinor(vcfInfo)
            indelCheck, type = checkIndel(position, gap_s, minority_depth, min_s, minorityVariant)
        if indelCheck == False: #If currect position is not indel, print previous
            if prev_position != None:
                newVCFlist.append(prev_position)
            if indellist != []: #Indel has ended, new position
                if prev_type == "insertion":
                    indelsequence = ""
                    nucleotideposition = indellist[0][1]
                    for indelposition in indellist:
                        indelsequence += indelposition[4]
                    for i in range(len(indellist)):
                        if i == 0: #Add sequence
                            indellist[i][4] = indelsequence
                            indellist[i][3] = indelsequence[0]
                        else:
                            indellist[i][4] = indelsequence[0:-i]
                            indellist[i][3] = indelsequence[0]
                        indellist[i][1] = nucleotideposition
                        indelvcfInfo = indellist[i][7].split(";") #add DP2
                        ad6 = indelvcfInfo[5][4:]
                        ad6list = ad6.split(",")
                        ad6listint = []
                        for number in ad6list:
                            ad6listint.append(int(number))
                        dp2 = max(ad6listint)
                        indelvcfInfo[-1] = "DP2=" + str(dp2)
                        indelvcfInfo = ";".join(indelvcfInfo)
                        indellist[i][7] = indelvcfInfo
                        indellist[i][5] = "0"
                    for item in indellist:
                        newVCFlist.append(item)

                elif prev_type == "deletion":
                    indelsequence = ""
                    nucleotideposition = indellist[0][1]
                    for indelposition in indellist:
                        indelsequence += indelposition[3]
                    for i in range(len(indellist)):
                        if i == 0:  # Add sequence
                            indellist[i][3] = indelsequence
                            indellist[i][4] = indelsequence[0]
                        else:
                            indellist[i][3] = indelsequence[0:-i]
                            indellist[i][4] = indelsequence[0]
                        indellist[i][1] = nucleotideposition
                        indelvcfInfo = indellist[i][7].split(";")  # add DP2
                        ad6 = indelvcfInfo[5][4:]
                        ad6list = ad6.split(",")
                        deletions = int(ad6list[-1])
                        dp = int(indelvcfInfo[0][3:])
                        indelvcfInfo[-1] = "DP2=" + str(deletions)
                        indelvcfInfo[1] = "AD=" + str(deletions)
                        indelvcfInfo[2] = "AF=" + str(deletions/dp)[0:4]
                        indelvcfInfo = ";".join(indelvcfInfo)
                        indellist[i][7] = indelvcfInfo
                        indellist[i][5] = "0"
                    for item in indellist:
                        newVCFlist.append(item)
                indellist = [] #Reset indellist

        elif indelCheck == True:
            if indellist == []: #First indel extension found
                indellist.append(prev_position)
                indellist.append(position)
            else:
                indellist.append(position)
        prev_type = type
        prev_position = position
    newVCFlist.append(prev_position)
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
    minority_count = int(sort_variantCount[-2])
    return minorityVariant, minority_depth, minority_count

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
    return positionType

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

