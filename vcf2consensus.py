#Import Libraries

import sys
import os
import argparse
import gzip



parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-ref', action="store", type=str, required=True, dest='referencefile', default="", help='Reference file')
parser.add_argument('-vcf', action="store", type=str, required=True, dest='vcf', default="", help='vcf')
parser.add_argument('-maj_s', action="store", type=float, dest='maj_s', default=0.7, help='Support for accepting majority variant')
parser.add_argument('-d', action="store", type=int, dest='depth', default=50, help='Depth threshold for including a position')
parser.add_argument('-min_s', action="store", type=float, dest='min_s', default=0.15, help='Support for accepting minority variant')
parser.add_argument('-variantCode', action="store", type=int, dest='variantCode', default=1, help='1 = Majority positions only, 2 = Minority positions only, 3 = UPAC minority positions')

args = parser.parse_args()

referencefile = args.referencefile
vcf = args.vcf
variantCode = args.variantCode
depth = args.depth
maj_s = args.maj_s
min_s = args.min_s

#['A','C','G','T','N','-']
UPACdict = dict()
UPACdict["AG"] = "R"
UPACdict["CT"] = "Y"
UPACdict["CG"] = "S"
UPACdict["AT"] = "W"
UPACdict["GT"] = "K"
UPACdict["AG"] = "M"

 # SHOULD BE MADE TO HANDLE MULTIPLE BASES PER VCF POSITION!!

def main():
    checkInput(variantCode)
    header, sequence = loadSequence(referencefile)
    vcfList = loadVCF(vcf)
    header, sequence, insertDict = consensusMaker(header, sequence, vcfList, depth, maj_s, min_s, UPACdict)
    header, sequence = handleInserts(header, sequence, insertDict)
    print (header)
    print (sequence)

def handleInserts(header, sequence, insertDict):
    insertsequence = ""
    for insert in insertDict:
        for value in insertDict[insert]:
            insertsequence = insertsequence + value[4]
        sequence = sequence[0:int(insert)] + insertsequence + sequence[int(insert):]
        insertsequence = ""
    return header, sequence

def consensusMaker(header, sequence, vcfList, depth, maj_s, min_s, UPACdict):
    insertDict = dict()
    for position in vcfList:
        vcfInfo = position[7].split(";")
        if variantCode == 1:  #Majority positions
            if float(vcfInfo[0][3:]) >= depth and float(vcfInfo[2][3:]) >= maj_s: #CheckDepth and maj_s and do not correct insertions or deletions
                if position[1] != "0":
                    previousPositions = position
                    if position[4] == "<->": #Handle deletion
                        sequence[int(position[1])] = "-"
                    else:
                        sequence[int(position[1])] = position[4]
                else:
                    if position[3] != "<->" or position[4] != "<->":  #Do not include gaps already in the reference
                        if previousPositions[1] in insertDict:
                            insertDict[previousPositions[1]].append(position)
                        else:
                            insertDict[previousPositions[1]] = [position]
        if variantCode == 2:
            if float(vcfInfo[0][3:]) >= depth: #CheckDepth and maj_s and do not correct insertions or deletions
                variants = vcfInfo[5][4:]
                variantCount = variants.split(",")
                sort_variantCount = []
                for i in range(len(variantCount)):
                    sort_variantCount.append(int(variantCount[i]))
                sort_variantCount.sort()
                minority_index = variantCount.index(str(sort_variantCount[-2]))
                dp = int(vcfInfo[0][3:])
                variantList = ['A','C','G','T','N','-']
                minorityVariant = variantList[minority_index]
                minority_depth = int(sort_variantCount[-2]) / dp
                if minority_depth > min_s:
                    if position[1] != "0":
                        previousPositions = position
                        if position[4] == "<->": #Handle deletion
                            sequence[int(position[1])] = "-"
                        else:
                            sequence[int(position[1])] = minorityVariant
                    else:
                        if position[3] != "<->" or position[4] != "<->":  #Do not include gaps already in the reference
                            if previousPositions[1] in insertDict:
                                insertDict[previousPositions[1]].append(position)
                            else:
                                insertDict[previousPositions[1]] = [position]
        if variantCode == 3:
            if float(vcfInfo[0][3:]) >= depth:  # CheckDepth and maj_s and do not correct insertions or deletions
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
                majorityVariant = position[4].upper()
                if minority_depth > min_s:
                    if position[1] != "0":
                        previousPositions = position
                        if position[4] == "<->":  # Handle deletion
                            sequence[int(position[1])] = "-"
                        else:
                            if (majorityVariant + minorityVariant) in UPACdict:
                                sequence[int(position[1])] = UPACdict[majorityVariant + minorityVariant]
                            elif (minorityVariant + majorityVariant) in UPACdict:
                                sequence[int(position[1])] = UPACdict[minorityVariant + majorityVariant]
                    else:
                        if position[3] != "<->" or position[4] != "<->":  # Do not include gaps already in the reference
                            if previousPositions[1] in insertDict:
                                insertDict[previousPositions[1]].append(position)
                            else:
                                insertDict[previousPositions[1]] = [position]
    sequence = ("").join(sequence)
    return header, sequence, insertDict

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
    return vcfList

def checkInput(variantCode):
    if variantCode > 3 or variantCode < 1:
        sys.exit("VariantCode must be either 1, 2 or 3")



if __name__ == '__main__':
    main()