import os
import sys
import argparse

from .bin_scanner import Scanner
from .bedsig import SimSeq


AP = argparse.ArgumentParser(
        description="rRNA depletion ratio evaluator",
        formatter_class=argparse.RawTextHelpFormatter,
        )

AP_subparsers = AP.add_subparsers(
        help="Sub-commands (use with -h for more info)"
        )

######################################################################################
##### Depletion ratio
######################################################################################
def _dep_rate(args):
    """Calculate depletion ratio of each bin."""
    scanner = Scanner(args)
    scanner.rbatch()

Dep_rate = AP_subparsers.add_parser('deprate', help=_dep_rate.__doc__)
Dep_rate.add_argument('-cfg', metavar='config', help='Config file, refer to example for details', required=True)
Dep_rate.add_argument('-p', metavar='process', help='Number of processor, default=4', default=4, type=int)
Dep_rate.set_defaults(func=_dep_rate)


######################################################################################
##### Simulating reads
######################################################################################
def _mock_seqs(args):
    """Simulating sequence reads from a reference genome and bed file.
     It can simulate diploid genomes with SNP and indels based on specific
     database (e.g., 1000 genomes)"""
    mocker = SimSeq(args)
    mocker.worker()

Mock = AP_subparsers.add_parser('mock', help=_mock_seqs.__doc__)
Mock.add_argument('-cfg', metavar='config', help='Config file, refer to example for details', required=True)
Mock.add_argument('-p', metavar='process', help='Number of processor, default=4', default=1, type=int)
Mock.set_defaults(func=_mock_seqs)


def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)
