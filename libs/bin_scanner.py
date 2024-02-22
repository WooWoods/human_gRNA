import os
import sys
import argparse
from multiprocessing import Pool

from .utils import read_config
from .process_bin import BinBuilder
from .process_bam import BamParser


class Scanner:
    def __init__(self, args):
        fcfg = args.cfg
        self.config = read_config(fcfg)
        self.p = args.p

        self.fh = open('rRNA_dep_rate2.tsv', 'w')
        
        bam_dir = self.config['bam_dir']

        self.bam_ref = os.path.join(bam_dir, self.config['bam_ref'])
        self.bam_del = os.path.join(bam_dir, self.config['bam_del'])

        ref_rRNA_ratio = self.config['ref_rRNA_ratio']
        del_rRNA_ratio = self.config['del_rRNA_ratio']

        self.ref_ratio = 1 - float(ref_rRNA_ratio)
        self.del_ratio = 1 - float(del_rRNA_ratio)

    def rbatch(self):
        bam_con = BamParser(self.bam_ref)
        bam_del = BamParser(self.bam_del)

        self.ref_reads = self.config['ref_reads'] or bam_con.total_reads()
        self.del_reads = self.config['del_reads'] or bam_del.total_reads()
        bam_con.close()
        bam_del.close()

        bin_obj = BinBuilder(self.config['rRNA_interval'])
        regions = bin_obj.interval_batch(self.p)

        pool = Pool(processes=self.p)
        results = []
        for region in regions:
            results.append(pool.apply_async(func=scanner, args=(self.bam_ref, self.bam_del, region,)))
        pool.close()
        pool.join()

        for res in results:
            self.recorder(res.get())

        self.fh.close()

    def recorder(self, results):
        for res in results:
            arr = res
            *region ,n_con, n_del = res
            fpkm_con = fpkm(n_con, int(self.ref_reads))
            fpkm_del = fpkm(n_del, int(self.del_reads))
            fpkm_con_normed = self.normalize(fpkm_con)

            ratio = del_ratio(fpkm_con_normed, fpkm_del)

            arr.extend([fpkm_con_normed, fpkm_del, ratio])

            arr_out = list(map(str, arr))
            print(arr_out)
            self.fh.write('\t'.join(arr_out) + '\n')

    def normalize(self, fpkm_con):
        return fpkm_con * self.del_ratio / self.ref_ratio

def scanner(bam_con, bam_del, region):
    results =[]

    bam_con_parser = BamParser(bam_con)
    bam_del_parser = BamParser(bam_del)

    for bin_ in region:
        n_con = bam_con_parser.bin_count(bin_)
        n_del = bam_del_parser.bin_count(bin_)
        print(bin_ + [n_con, n_del])
        results.append(bin_ + [n_con, n_del])
    bam_con_parser.close()
    bam_del_parser.close()
    return results

def del_ratio(fpkm_con, fpkm_del):
    return (fpkm_con - fpkm_del)/fpkm_con

def fpkm(reads, total_reads, bin_len=100):
    return reads * 1e9/(bin_len * total_reads)
    

