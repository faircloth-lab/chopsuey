#!/usr/bin/env python
# encoding: utf-8
"""
File: chopsuey.py
Author: Brant Faircloth

Created by Brant Faircloth on 06 August 2012 14:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import pdb
import glob
import string
import argparse
from Bio import SeqIO
from Bio import Restriction

import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "dna",
            action=FullPaths,
            help="""A folder containing your DNA source(s) of interest"""
        )
    parser.add_argument(
            "enzymes",
            type=str,
            args='+',
            help="""The enzyme(s) to chop with""",
        )
    return parser.parse_args()

def get_all_fasta_files(dna):
    fastas = []
    for ext in ['*.fasta', '*.fsa', '*.fa']:
        fastas.extend(glob.glob(os.path.join(dna, ext)))
    return fastas

def main():
    args = get_args()
    fastas = get_all_fasta_files(args.dna)
    for source in fastas:
        seq = SeqIO.read(open(file,'rU'),'fasta')
        for attr, value in Restriction.__dict__.iteritems():
            pdb.set_trace()

if __name__ == '__main__':
    main()
