#!/usr/bin/env python3
#

"""
This script takes gene count table: gene_name_1 \t read_count
1) combines read counts if exons of one gene annotated separately
2) assigns gene sizes from another file
"""

import argparse
import re
from collections import defaultdict
import sys


__author__ = "Ekaterina Osipova, 2021."


def read_anno_bed(anno_file):
    ## Reads annotation file into a dict: transcript_name : transcript_size

    transcript_dict = defaultdict(list)
    with open(anno_file, 'r') as inf:
        for bed_line in inf.readlines():
            name = bed_line.split()[3]
            if len(bed_line.split()) == 12:
                block_sizes = [int(t) for t in bed_line.split()[10].rstrip(',').split(',')]
                transcript_dict[name] = block_sizes
            elif len(bed_line.split()) == 6 or len(bed_line.split()) == 4:
                start = int(bed_line.split()[1])
                stop = int(bed_line.split()[2])
                transcript_dict[name].append(stop - start)
            else:
                print('Provide a valid bed file!')
                sys.exit(1)
    return transcript_dict


def read_count_table(count_file):
    ## Reads read count file into a dict. Combines _2 exons in one transcript

    count_dict = defaultdict(list)
    with open(count_file, 'r') as inf:
        for line in inf.readlines():
            exon_name = line.split()[0]
            exon_count = int(line.split()[1])
            exon_suffix = r'_[0-9]+$'
            transc_name = re.sub(exon_suffix, '', exon_name)
            count_dict[transc_name].append(exon_count)
    return count_dict


def output_updated_count_table(count_dict, transcript_dict):
    ## Prints new count table: gene_name \t read_count \t gene_size

    for transc in count_dict:
        reads = sum(count_dict[transc])
        transc_size = transcript_dict[transc]
        if len(transc_size) > 0:
            count_line = '{}\t{}\t{}'.format(transc, reads, sum(transc_size))
        else:
            count_line = '{}\t{}\t{}'.format(transc, reads, 1)
        print(count_line)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--readcount', type=str, help='read count table: gene_name\tread_count')
    parser.add_argument('-a', '--anno', type=str, help='bed12 annotation file to extract gene sizes from')
    args = parser.parse_args()

    ## Read annotation file
    transcript_dict = read_anno_bed(args.anno)

    ## Read read count table
    count_dict = read_count_table(args.readcount)
    # print(count_dict)
    # sys.exit()

    ## Print updated count table
    output_updated_count_table(count_dict, transcript_dict)


if __name__ == "__main__":
    main()