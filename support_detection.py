#Import Libraries

import sys
import os
import argparse
import gzip
import json

#Identifies variants according to inputput type, and can give a binary output when searching for a specific variant

parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-p', action="store", type=str, required=True, nargs='+', dest='positions', default="", help='Positions')
parser.add_argument('-vcf', action="store", type=str, required=True, nargs='+', dest='vcf', default="", help='vcf. ')

args = parser.parse_args()

positions = args.positions
vcf = args.vcf

def main():
    position_list = []
    del_flag = False
    for file in vcf:
        vcfList = loadVCF(file)
        for item in vcfList:
            if item[1] not in positions and item[1] != "0":
                del_flag = False
            if item[1] in positions:
                del_flag = True
                position_list.append(item)
            elif del_flag == True and (item[3] != "<->" and item[4] != "<->"):
                position_list.append(item)
    support_object = dict()
    for position in position_list:
        vcfInfo = position[7].split(";")
        dp = int(vcfInfo[0][3:])
        ad6 = vcfInfo[5][4:]
        ad6list = ad6.split(",")
        gaps_d = int(ad6list[-1])
        if position[1] in support_object:
            support_object[position[1]].append(str(dp) + "-" + str(gaps_d))
        else:
            support_object[position[1]] = [str(dp) + "-" + str(gaps_d)]
    for position in support_object:
        total_d = []
        total_g = []
        support_list = []
        for i in range(len(support_object[position])):
            count = support_object[position][i].split("-")
            total_d.append(int(count[0]))
            total_g.append(int(count[1]))
            support_list.append(int(count[1])/int(count[0]))
        avg_s = sum(total_g)/sum(total_d)
        print ("{} has an average of {} support with a min of {} and max of {}".format(position, avg_s, min(support_list), max(support_list)))

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