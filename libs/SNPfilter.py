import os
import re
import sys

from bisect import bisect_left, insort_left
from collections import defaultdict



def load_interval(fbed):
    intervals = defaultdict(list)
    with open(fbed) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            arr = line.split()
            chrom = arr[0]
            start = int(arr[1])
            end = int(arr[2])
            insort_left(intervals[chrom], [start, end], key=lambda r: r[0])
    return intervals

def in_bed(chrom, pos, intervals):
    regions = intervals[chrom]
    index = bisect_left(regions, [pos, pos])
    for i in [index-1, index] :
        try:
            region = regions[i]
        except IndexError:
            continue
        if pos > region[0] and pos < region[1]:
            return True
    return False

def ambiguous(alt):
    m = re.search(r'(\d+)([ATCG]+)', alt)
    if m is None:
        return False
    if int(m.group(1)) > 1:
        return True

def filter_snp(SNPdb, intervals, cutoff=0.1):
    with open(SNPdb) as fh, open('SNP_candidates.tsv', 'w') as fo:
        for line in fh:
            arr = line.split()
            if float(arr[4]) < cutoff or float(arr[4]) > 1 - cutoff:
                continue
            if ambiguous(arr[3]):
                continue
            chrom = arr[0]
            pos = int(arr[1])
            if in_bed(chrom, pos, intervals):
                fo.write(line)


if __name__ == '__main__':
    fbed, SNPdb = sys.argv[1:]
    intervals = load_interval(fbed)
    filter_snp(SNPdb, intervals)





