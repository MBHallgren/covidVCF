#Import Libraries

import sys
import os
import argparse
import gzip

#Identifies variants according to inputput type, and can give a binary output when searching for a specific variant

parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-p', action="store", type=str, required=True, nargs='+', dest='positions', default="", help='Positions')
parser.add_argument('-vcf', action="store", type=str, required=True, dest='vcf', default="", help='vcf')
parser.add_argument('-maj_s', action="store", type=float, dest='maj_s', default=0.7, help='Support for accepting majority variant')
parser.add_argument('-min_s', action="store", type=float, dest='min_s', default=0.15, help='Support for accepting minority variant')
parser.add_argument('-d', action="store", type=int, dest='depth', default=50, help='depth threshold for including a position')
parser.add_argument('-f', action="store", type=str, dest='filter', default='1', help='Filter for either all variants (1), Majority only (2), Minority only(3)')
parser.add_argument('-o', action="store", type=str, dest='outputformat', default='1', help='Output format: Check all variants an given positions, and output those found (1), All positions must be found (2) ')


args = parser.parse_args()

positions = args.positions
vcf = args.vcf
depth = args.depth
maj_s = args.maj_s
min_s = args.min_s
filter = args.filter
outputformat = args.outputformat

def main():
    vcfList = loadVCF(vcf)
    checkForVariants(vcfList, positions, depth, maj_s, min_s, filter, outputformat)



def checkForVariants(vcfList, positions, depth, maj_s, min_s, filter, outputformat):
    found_positions = []
    for position in positions:
        majority = ""
        minority = ""
        for i in range(len(vcfList)):
            if position == vcfList[i][1]:
                vcfInfo = vcfList[i][7].split(";")
                variants = vcfInfo[5][4:]
                if int(vcfInfo[0][3:]) >= depth: #Check depth
                    if filter == '2': #Check support for majority:
                        if int(vcfInfo[1][3:])/int(vcfInfo[0][3:]) >= maj_s:
                            majority = "Position: {} <{}> is a Majority Variant with support: {}".format(position, vcfList[i][3], int(vcfInfo[1][3:])/int(vcfInfo[0][3:]))
                    if filter == '3': #Check support for majority:
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
                        if minority_depth >= min_s: #Support
                            minority = "Position: {} <{}> is a Minority Variant with support: {}".format(position,minorityVariant, minority_depth)
                    if filter == '1':
                        if int(vcfInfo[1][3:])/int(vcfInfo[0][3:]) >= maj_s:
                            majority = "Position: {} <{}> is a Majority Variant with support: {}".format(position, vcfList[i][3], int(vcfInfo[1][3:])/int(vcfInfo[0][3:]))
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
                        if minority_depth >= min_s:  # Support
                            minority = "Position: {} <{}> is a Minority Variant with support: {}".format(position,minorityVariant,minority_depth)
        if outputformat == '1':
            if majority != "":
                print (majority)
            if minority != "":
                print (minority)
        if outputformat == '2':
            if minority != "":
                found_positions.append("minority,{},{}".format(minorityVariant,position))
            elif majority != "":
                found_positions.append("majority,{},{}".format(vcfList[i][3], position))
    if outputformat == '2':
        if len(positions) == len(found_positions):
            print ("All positions were identified as variants")
            print ("{}/{} positions found".format(len(found_positions), len(positions)))
            print (found_positions)
        else:
            print("All positions were NOT identified as variants")
            print("{}/{} positions found".format(len(found_positions), len(positions)))
            print(found_positions)










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
    infile.close()
    return vcfList

if __name__ == '__main__':
    main()