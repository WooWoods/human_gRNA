import os
import sys
import random

from bisect import bisect_left, insort_left
from collections import defaultdict



class BinBuilder:
    def __init__(self, bed, bin_len= 100, overlap=50):
        self.bin_len = bin_len
        self.overlap = overlap
        self.bed = bed
    
    def interval(self):
        allbins = []
        try:
            fh = open(self.bed)
        except Exception as err:
            print(f"Open file error: {err}")
            raise

        for line in fh:
            if line.startswith('@'):
                continue
            arr = line.split()
            bins = self.splitter(arr[0:3])
            allbins.extend(bins)
        fh.close()
        return allbins

    def splitter(self, region):
        bins = []
        chrom = region[0]
        start = int(region[1])
        end = int(region[2])

        bin_start = start
        while bin_start < end:
            bin_end = bin_start + self.bin_len
            bins.append([chrom, bin_start, bin_end])
            bin_start += self.overlap

        return bins

    def interval_batch(self, p):
        intervals = []
        allbins = self.interval()

        n_bins = len(allbins)
        batch_lines = round(n_bins/p)
    
        pos = 0
        for i in range(p):
            if i == p - 1:
                intervals.append(allbins[pos:])
            else:
                intervals.append(allbins[pos: pos+batch_lines])
            pos += batch_lines
        return intervals


def random_bins(fbed, bin_len=100):
    intervals = defaultdict(list)
    steps = bin_len * 0.4
    with open(fbed) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            arr = line.split()
            chrom = arr[0]
            start = int(arr[1]) - 1
            end = int(arr[2])
            rlen = end - start
            if rlen < bin_len:
                insort_left(intervals[chrom], [start, end], key=lambda r: r[0])
            else:
                segs = round(rlen/steps)
                for i in range(segs):
                    pos = random.randint(start, end)
                    pos_end = min(end, pos+bin_len)
                    if pos_end - pos < 30:
                        continue
                    insort_left(intervals[chrom], [pos, pos_end], key=lambda r: r[0])
    return intervals




