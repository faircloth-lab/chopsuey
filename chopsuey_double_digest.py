#!/usr/bin/env python
# encoding: utf-8
"""
File: chopsuey_double_digest.py
Author: Brant Faircloth

Created by Brant Faircloth on 06 August 2012 14:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Return summary data from a double digest of a directory of
targets using one or more restriction enzymes.

Usage:

    python chopsuey_double_digest.py /folder/with/some/fastas PstI BfaI NdeI AseI MseI

This will return combinations of the first enzyme (PstI) with the other enzymes.

"""

import os
import pdb
import glob
import math
import numpy
import argparse
from Bio import SeqIO
from Bio import Restriction

#import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Summary data from a double digest of a directory of
targets""")
    parser.add_argument(
            "dna",
            action=FullPaths,
            help="""A folder containing your DNA source(s) of interest"""
        )
    parser.add_argument(
            "first",
            type=str,
            help="""The first enzyme you will cut with""",
        )
    parser.add_argument(
            "others",
            type=str,
            nargs='+',
            help="""The additional enzyme(s) you will cut with""",
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
    print """
  ___  _  _   __  ____    ____  _  _  ____  _  _
 / __)/ )( \ /  \(  _ \  / ___)/ )( \(  __)( \/ )
( (__ ) __ ((  O )) __/  \___ \) \/ ( ) _)  )  /
 \___)\_)(_/ \__/(__)    (____/\____/(____)(__/

         Virtual Enzyme Double Digest
"""
    for fasta in fastas:
        fasta_name = os.path.basename(fasta)
        print "\nWorking on {}...".format(fasta_name)
        fasta_results = {}
        seq = SeqIO.read(open(fasta, 'rU'), 'fasta')
        print "\tCutting\n\t======================"
        for attr, value in Restriction.__dict__.iteritems():
            if (attr in args.first) or (attr in args.others):
                cut_site = value.elucidate()
                cuts = sorted(value.search(seq.seq))
                print "\tEnzyme {}...".format(attr)
                fasta_results[attr] = {
                        'cut_site': cut_site,
                        'cuts': cuts,
                        'seq_len': len(seq.seq)
                    }
        cut_results = {}
        base = fasta_results[args.first]
        print "\tDouble-cutting\n\t======================"
        for enzyme, values in fasta_results.iteritems():
            all_cuts = base['cuts']
            if enzyme != args.first:
                all_cuts.extend(values['cuts'])
                all_cuts = sorted(all_cuts)
                lengths = [v - all_cuts[k] for k, v in enumerate(all_cuts[1:])]
                #pdb.set_trace()
                combo = "{0}-{1}".format(args.first, enzyme)
                print "\tCombo {}...".format(combo)
                cut_results[combo] = {
                        'cuts': len(all_cuts),
                        'cuts_per_bp': float(len(all_cuts)) / values['seq_len'],
                        'cut_length': numpy.mean(lengths),
                        'cut_95_ci': 1.96 * numpy.std(lengths, ddof=1) / math.sqrt(len(lengths)),
                        'min': min([i for i in lengths if i > 0]),
                        'max': max(lengths)
                    }
        results[fasta_name] = cut_results
        #pdb.set_trace()
    print "\nFile,Enzyme Combo,Cuts,Cuts per bp,Avg. fragment length,95CI cut length,Min(size),Max(size)"
    for fasta_name, values in results.iteritems():
        for combo_name, data in values.iteritems():
            print "{0},{1},{2[cuts]},{2[cuts_per_bp]},{2[cut_length]},{2[cut_95_ci]},{2[min]},{2[max]}".format(
                    fasta_name,
                    combo_name,
                    data)


if __name__ == '__main__':
    main()
