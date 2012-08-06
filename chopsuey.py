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
            nargs='+',
            help="""The enzyme(s) to chop with""",
        )
    return parser.parse_args()

def get_all_fasta_files(dna):
    fastas = []
    for ext in ['*.fasta', '*.fsa', '*.fa', '*.fas']:
        fastas.extend(glob.glob(os.path.join(dna, ext)))
    return fastas

def main():
    args = get_args()
    fastas = get_all_fasta_files(args.dna)
    results = {}
    for fasta in fastas:
        fasta_name = os.path.basename(fasta)
        fasta_results = {}
        seq = SeqIO.read(open(fasta, 'rU'), 'fasta')
        for attr, value in Restriction.__dict__.iteritems():
            if attr in args.enzymes:
                name = attr
                cut_site = value.elucidate()
                cuts = sorted(value.search(seq.seq))
                lengths = [v - cuts[k - 1] for k, v in enumerate(cuts[1:])]
                fasta_results[name] = {
                        'cut_site': cut_site,
                        'cuts': len(cuts),
                        'cuts_per_bp': float(len(cuts)) / len(seq.seq),
                        'cut_length': sum(lengths) / len(lengths)
                    }
        results[fasta_name] = fasta_results
    for fasta_name, values in results.iteritems():
        for enzyme_name, data in values.iteritems():
            print "{0},{1},{2[cuts]},{2[cuts_per_bp]},{2[cut_site]},{2[cut_length]}".format(
                    fasta_name,
                    enzyme_name,
                    data)
    pdb.set_trace()

if __name__ == '__main__':
    main()
