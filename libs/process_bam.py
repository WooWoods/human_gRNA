import os
import re
import sys

import pysam

def process_CIGAR(pos, cigar):
    block_lengths  = [] 
    soft_clippings = []
    blocks = 0
    for (ctype, length) in cigar:
        if ctype == 0: #M
            blocks += length
        elif ctype == 2: #D
            blocks += length
        elif ctype == 3: #N
            block_lengths.append((pos, pos + blocks))
            pos += (blocks + length)
            blocks = 0
        elif ctype == 4: #S
            soft_clippings.append(length)
    block_lengths.append((pos, pos + blocks))
    return block_lengths, soft_clippings

class BamParser:
    def __init__(self, bam, maxSoftClip=50):
        self.maxSoftClip = maxSoftClip
        self.bam = bam
        self.bamfile = pysam.AlignmentFile(bam, 'rb')

    def total_reads(self):
        pattern = re.compile(r'(\d+).*in total')
        res = pysam.flagstat(self.bam)
        m = pattern.search(res)
        total_reads = int(m.group(1))
        return total_reads
        
    def bin_count(self, region):
        n = 0
        chrom, start, end = region
        for read in self.bamfile.fetch(chrom, start, end):
            r_blocks, soft_clippings = process_CIGAR(read.pos, read.cigar)

            if any(SC > self.maxSoftClip for SC in soft_clippings):
                continue
            
            if self.valid_read(region, r_blocks):
                n += 1
        return n

    def valid_read(self, region, r_blocks):
        _, start, end = region
        bin_len = end - start
        m_len = 0
        for bs, be in r_blocks:
            if bs > end or be < start:
                continue
            if be < end:
                st = max(bs, start)
                m_len = be - st
            if bs > start:
                ed = min(end, be)
                m_len = ed - bs

            if m_len/bin_len > 0.5:
                return True
        return False

    def close(self):
        self.bamfile.close()

