#Import Libraries

import sys
import os
import argparse
import gzip

#This script is for filtering variants from a VCF file. Specify a set of positions, and the type of variants to be filtered out.

parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-p', action="store", type=str, required=True, nargs='+', dest='positions', default="", help='Positions')
parser.add_argument('-vcf', action="store", type=str, required=True, nargs='+', dest='vcf', default="", help='vcf(s)')
parser.add_argument('-vtype', action="store", type=str, required=True, dest='vtype', default="a", help='[all, min, major] = [all position variants, minority variants in given positions, majority variants in given positions')
parser.add_argument('-min_threshold', action="store", type=float, dest='min_threshold', default=0.3, help='Minority call threshold')
parser.add_argument('-maj_threshold', action="store", type=float, dest='maj_threshold', default=0.7, help='Majority call threshold')
parser.add_argument('-d', action="store", type=int, dest='depth', default=50, help='depth threshold for including a position')

args = parser.parse_args()

vtype = args.vtype
positions = args.positions
vcf = args.vcf
depth = args.depth
if vtype == "min":
    threshold = args.min_threshold
elif vtype == "maj":
    threshold = args.maj_threshold
else:
    threshold = 0.5


def main():
    callVariants(vcf, positions, vtype, threshold, depth)

def minorityCall(list, threshold, depth):
    info = list[7]
    infolist = info.split(";")
    af = float(infolist[2][3:])
    ad = float(infolist[0][3:])
    if threshold > af and ad > depth: #Minimum allele freq allowed for minorities
        print (af)
        return True
    else:
        return False

def majorityCall(list, threshold, depth):
    info = list[7]
    infolist = info.split(";")
    af = float(infolist[2][3:])
    ad = float(infolist[0][3:])
    if af > threshold and ad > depth: #maximum allele freq allowed for majorities
        print (af)
        return True
    else:
        return False

def callVariantType(vtype, line, threshold, depth):
    if vtype == "all":
        list = line.split("\t")
        info = list[7]
        infolist = info.split(";")
        ad = float(infolist[0][3:])
        if ad > depth:
            print (line)
            returnValue = True
    elif vtype == "min":
        returnValue = minorityCall(line.split("\t"), threshold, depth)
        if returnValue:
            print (line)
    elif vtype == "maj":
        returnValue = majorityCall(line.split("\t"), threshold, depth)
        if returnValue:
            print (line)
    else:
        sys.exit("vtype (variantType) was not corrently given [all, min, maj]")
    return returnValue


def callVariants(vcf, positions, vtype, threshold, depth):
    insertflag = False
    for file in vcf:
        if file[-2:] == 'gz':
            infile = gzip.open(file, 'rb')
            for line in infile:
                line = line.decode()
                line = line.rstrip()
                if line[0] == '#':
                    print(line)
                else:
                    line = line.rstrip()
                    linelist = line.split("\t")
                    if insertflag == True:
                        if linelist[1] == "0":
                            print(line)
                        else:
                            insertflag = False
                    if linelist[1] in positions:
                        returnCall = callVariantType(vtype, line, threshold, depth)
                        if returnCall:
                            insertflag = True

        else:
            infile = open(file, 'r')
            for line in infile:
                line = line.rstrip()
                if line[0] == '#':
                    print(line)
                else:
                    line = line.rstrip()
                    linelist = line.split("\t")
                    if insertflag == True:
                        if linelist[1] == "0":
                            print(line)
                        else:
                            insertflag = False
                    if linelist[1] in positions:
                        returnCall = callVariantType(vtype, line, threshold, depth)
                        if returnCall:
                            insertflag = True
        infile.close()


if __name__ == '__main__':
    main()

