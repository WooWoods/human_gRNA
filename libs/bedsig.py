import os
import re
import sys

from uuid import uuid4
from random import random, randint, choice
from collections import defaultdict
from bisect import bisect_right, insort_left

from .utils import read_config
from .process_bin import random_bins




def fa_parser(fa):
    header = ''
    seq = []
    with open(fa) as fh:
        for line in fh:
            if line.startswith('>'):
                if seq:
                    yield (seqid, ''.join(seq))
                header = line.strip()[1:]
                seqid, *tmp = header.split()
                seqid = re.sub('chr', '', seqid)
                seq = []
            else:
                seq.append(line.strip())
        yield (seqid, ''.join(seq))

class SimSeq:
    def __init__(self, args):
        fcfg = args.cfg
        self.config = read_config(fcfg)
        self.p = args.p
        self.fa = self.config['genome']
        self.n = int(self.config['n'])
        self.snps = self.config['SNPdb']
        self.add_snp = self.config['addSNP']
        fbed = self.config['rRNA_interval']
        self.random_intervals = random_bins(fbed)
        self.load_snps()

    def load_snps(self):
        self.snp_db = defaultdict(list)
        self.search_db = defaultdict(list)
        with open(self.snps) as fh:
            for line in fh:
                arr = line.split()
                chrom = arr[0]
                start = int(arr[1])

                key = tuple(arr[0:2])
                self.snp_db[key].append(arr[2:5])
                insort_left(self.search_db[chrom], start)
        #print(self.snp_db)

    def snps_in_region(self, chrom, region):
        res = []
        db = self.search_db[chrom]
        first = bisect_right(db, region[0])
        last = bisect_right(db, region[1])

        idx1 = max(0, first-1)
        for pos in db[idx1:last]:
            if  pos > region[0] and pos < region[1]:
                res.append((chrom, str(pos)))
        return res

    def seq_with_var(self, read, region, pos, snp_info):
        pattern = re.compile(r'(\d+)([ATCG]+)')
        ref, alt, freq = snp_info
        pos = int(pos)

        m = pattern.search(alt)
        alt_pos = pos - region[0] - 1

        # SNP or deletion
        if m is None:
            try:
                # deletion
                tmp = int(alt)
                new_read = read[0:alt_pos] + read[alt_pos+len(ref):]
            except ValueError:
                new_read = read[0:alt_pos] + alt + read[alt_pos+1:]
        # multi alts
        else:
            alt = m.group(2)
            new_read = read[0:alt_pos] + alt + read[alt_pos+len(alt):]

        return new_read

    def worker(self):
        fh = open('mock.fa', 'w')
        for chrom, seq in fa_parser(self.fa):
            try:
                regions = self.random_intervals[chrom]
                if not regions:
                    continue
                for region in regions:
                    snps = self.snps_in_region(chrom, region)
                    read = seq[region[0]:region[1]]
                    #read_reverse = read_rc(read)
                    if self.add_snp == 'False':
                        for i in range(self.n):
                            seq_id = gen_seq_id(chrom, region)
                            fh.write(f'>{seq_id}\n{read}\n')
                            #fh.write(f'>{seq_id}_rc\n{read_reverse}\n')
                    else:
                        if not snps:
                            for i in range(self.n):
                                seq_id = gen_seq_id(chrom, region)
                                fh.write(f'>{seq_id}\n{read}\n')
                                #fh.write(f'>{seq_id}_rc\n{read_reverse}\n')
                        else:
                            snp = choice(snps)
                            _, pos = snp
                            snp_info = choice(self.snp_db[snp])
                            freq = float(snp_info[2])
                            output_read = self.seq_with_var(read, region, pos, snp_info)
                            #print(region)
                            #print(pos, snp_info)

                            for i in range(self.n):
                                prob = random()
                                if prob < freq:
                                    out_read = self.seq_with_var(read, region, pos, snp_info)
                                else:
                                    out_read = read

                                #read_reverse = read_rc(out_read)
                                seq_id = gen_seq_id(chrom, region)
                                fh.write(f'>{seq_id}\n{out_read}\n')
                                #fh.write(f'>{seq_id}_rc\n{read_reverse}\n')
            except KeyError:
                continue

        fh.close()

def gen_seq_id(chrom, region):
    num = randint(1000, 3000)
    string = str(uuid4())
    string_seg, *tmp = string.split('-')
    return f'{chrom}_{region[0]}_{region[1]}_{num}{string_seg}'

def read_rc(read):
    base_pairs = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'N': 'N'
            }
    return ''.join([base_pairs[b] for b in read][::-1])





