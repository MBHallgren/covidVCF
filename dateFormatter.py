#Import Libraries

import sys
import os
import argparse
import gzip
import json

#SNP_csv is output from variantIdentifier.py when it has been run with multiple VCF's at once
#date_csv should be a csv list with the filenames of the VCF and sampling date (file1,02-03-2021)

parser = argparse.ArgumentParser(description='vcfFilter')
parser.add_argument('-snp_csv', action="store", type=str, required=True, dest='snp_csv', default="", help='snp_csv')
parser.add_argument('-date_csv', action="store", type=str, required=True, dest='date_csv', default="", help='date_csv')

args = parser.parse_args()

def main():
    snp_csv_list = loadCSV(args.snp_csv)
    date_csv_list = loadCSV(args.date_csv)
    date_dict = dateDict(date_csv_list)

def insertDate(snp_csv_list, date_dict):
    newlist = []
    for item in snp_csv_list:
        if

def dateDict(date_csv_list):
    date_dict = dict()
    for item in date_csv_list:
        if item[0] not in date_dict:
            date_dict[item[0]] = item[1]
        else:
            pass
    return date_dict


def loadCSV(file):
    csvList = []
    infile = open(file, 'r')
    for line in infile:
        csvList.append(line.rstrip().split(","))
    return csvList

if __name__ == '__main__':
    main()
