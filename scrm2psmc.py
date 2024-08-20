#!/usr/bin/python
# -*-coding:Latin-1 -*

import bisect
import sys
import textwrap

from time import time

from collections import defaultdict
from math import floor, ceil
from argparse import ArgumentParser, FileType, RawTextHelpFormatter


def get_configuration():

    parser = ArgumentParser(description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('-s', '--size', metavar='M', type=int, default=100,
        help="set the window size\n"
        "(default: %(default)s)")

    parser.add_argument('-o', dest='outputfile', metavar='FILE', type=FileType('w'),
        help="Output file\n"
        "(default: console)")

    parser.add_argument('inputfile', metavar='FILE', type=FileType('r'), nargs=1,
        help = "Input file: output of scrm (or VCF if option -v is provided) to convert into the input for PSMC (use `-` for STDIN)")

    parser.add_argument('-v', '--vcf', action='store_true',
                        help="Convert a VCF file (with genotypes from a single individual) into the input for PSMC\n"
                              "(default: False)")

    parser.add_argument('-m', dest='inputmaskfile', metavar='FILE', type=FileType('r'), nargs=1, required=True,
        help = "Input BED file with genomic regions that pass QC")

    parser.add_argument('-g', dest='genomeLength', metavar='FILE', type=FileType('r'), nargs=1, required=True,
        help = "Tab separated file with the chromosome name and its length, each line corresponding to a different chromosome")

    opts = parser.parse_args()
    opts.inputfile = opts.inputfile[0]
    opts.inputmaskfile = opts.inputmaskfile[0]
    opts.genomeLength = opts.genomeLength[0]
    return opts


def read_genome_size(fh):
    chrom_len = defaultdict()

    for line in fh:
        chrom, length = line.rstrip().split('\t')
        chrom_len[chrom] = int(length)

    return chrom_len


def read_scrm(chrom_len, fh=sys.stdin):
    sites = defaultdict()

    start = defaultdict()
    end = defaultdict()

    # Assuming that there is only one simulation and that it should be split into mutliple chromosomes
    previous=0;
    for chrom in chrom_len:
        start[chrom] = previous + 1
        end[chrom] = start[chrom] + chrom_len[chrom] - 1
        previous=end[chrom]

    for line in fh:
        # Skip the header
        if not line.startswith('positions:'):
            continue

        raw_sites = line.rstrip().split()

        for i in range(1, len(raw_sites)):
            position = round(float(raw_sites[i]))
            for chrom in chrom_len:
                if position >= start[chrom] and position <= end[chrom]:
                    try:
                        sites[chrom].append(position-start[chrom]+1)
                    except KeyError:
                        sites[chrom] = [position-start[chrom]+1]
    return sites # coordinates are 1-based


def read_vcf(fh=sys.stdin):
    sites = defaultdict()

    first_pos = defaultdict()
    last_pos = defaultdict()

    for line in fh:
        # Skip the header
        if line.startswith('#'):
            continue

        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, GT, *opt = \
            line.rstrip().split('\t')

        if not CHROM in sites:
            first_pos[CHROM]=int(POS)

        if GT[0]!=GT[2]:
            try:       
                sites[CHROM].append(int(POS))
            except KeyError:
                sites[CHROM]=[int(POS)]

        last_pos[CHROM]=int(POS)

    return sites, first_pos, last_pos
        

def read_BED(fh):
    intervals = defaultdict()

    for line in fh:
        # Skip the header
        if line.startswith('#'):
            continue

        chrom, start, end, *opt = \
            line.rstrip().split('\t')

        try:
            intervals[chrom].append((int(start), int(end)))
        except KeyError:
            intervals[chrom]=[(int(start), int(end))]

    return intervals


def getOverlap(a, b):
#    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    return max(0, min(a[1], b[1]) - max(a[0], b[0]+1) + 1) # b[0]+1 because b is 0-based while a is 1-based (but not for b[1], because b[1]+1 is not included)


def remove_hets_sorted_list(lst, value):
    while len(lst) > 0 and lst[0] < value:
        lst.pop(0)
    return lst


def remove_intervals_sorted_list(intervals, value):
    while len(intervals) > 0 and intervals[0][1] < value:
        intervals.pop(0)
    return intervals


def filter_sorted_list_by_intervals(lst, intervals):
    result = []
    i = j = 0
    while i < len(lst) and j < len(intervals):
        value = lst[i]
        start, end = intervals[j] # These coordinates should be 0-based and intervals should be [start, end)
        if value-1 < start: # value-1 because these are 1-based
            i += 1
        elif value-1 >= end:
            j += 1
        else:
            result.append(value)
            i += 1
    return result


def main():

    options = get_configuration()

    chrom_length = read_genome_size(options.genomeLength)
    if options.vcf:
        hets, first_pos, last_pos = read_vcf(options.inputfile)
    else:
        hets = read_scrm(chrom_length, options.inputfile)
    masks = read_BED(options.inputmaskfile)
    window_size = options.size

    for chrom in hets:
        print('>{}'.format(chrom), file=options.outputfile)

        sequence = str()

        hets[chrom].sort()
        try:
            masks[chrom].sort()
        except KeyError:
            continue

        hets[chrom] = filter_sorted_list_by_intervals(hets[chrom], masks[chrom])

        if len(hets[chrom])==0:
            continue

        if options.vcf:
            start=first_pos[chrom]
            end=last_pos[chrom]
        else:
            start=min(hets[chrom])
            end=max(hets[chrom])

        pointer_het = 0
        pointer_mask = 0

        for i in range(floor(start/window_size)*window_size, ceil(end/window_size)*window_size, window_size):
            covered = 0

            while hets[chrom][pointer_het] <= i and pointer_het < len(hets[chrom])-1: #not including the first position which should be from the previous window (i.e., the first window does not start at 0, but at 1)
                pointer_het += 1
            while masks[chrom][pointer_mask][1] <= i and pointer_mask < len(masks[chrom])-1:
                pointer_mask += 1

            pointer_overlapping_mask = int(pointer_mask)
            while masks[chrom][pointer_overlapping_mask][0] < i + window_size and pointer_overlapping_mask < len(masks[chrom])-1:
                covered += getOverlap((i + 1, i + window_size), masks[chrom][pointer_overlapping_mask])
                pointer_overlapping_mask += 1

            if covered < 0.1*window_size:
                sequence += "N"
                continue

            if hets[chrom][pointer_het] <= i + window_size and hets[chrom][pointer_het] > i:
                sequence += "K"
            else:
                sequence += "T"

        sequence = sequence.strip("N")

        wrapped_sequence = textwrap.wrap(sequence, 60)
        for line in wrapped_sequence:
            print(line, file=options.outputfile)

if __name__ == '__main__':
    main()


