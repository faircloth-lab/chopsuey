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
#import pdb
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
    parser.add_argument(
            "--opposite-ends",
            dest='opposite_ends',
            action="store_true",
            default=False,
            help="""Ensure one end goes with one enzyme, the other end with the other""",
        )
    parser.add_argument(
            "--filter",
            type=int,
            default=None,
            nargs='+',
            help="""Filter the resulting fragment size distribution""",
        )
    return parser.parse_args()


def get_all_fasta_files(dna):
    fastas = []
    for ext in ['*.fasta', '*.fsa', '*.fa', '*.fas']:
        fastas.extend(glob.glob(os.path.join(dna, ext)))
    return fastas


def get_single_digest(arguments, seq):
    single_cut_results = {}
    print "\tCutting\n\t======================"
    for attr, value in Restriction.__dict__.iteritems():
        # work only with enzymes input by user
        if (attr in arguments.first) or (attr in arguments.others):
            cut_site = value.elucidate()
            cuts = sorted(value.search(seq.seq))
            print "\tEnzyme {}...".format(attr)
            single_cut_results[attr] = {
                    'cut_site': cut_site,
                    'cuts': cuts,
                    'seq_len': len(seq.seq)
                }
    return single_cut_results


def get_opposite_end_results(arguments, seq, base, enzyme, compare, all_cuts, cut_results):
    fragments = []
    base_sites = base['cut_site'].replace('_', '').split('^')
    compare_sites = compare['cut_site'].replace('_', '').split('^')
    ends = ([base_sites[1], compare_sites[0]], [compare_sites[1], base_sites[0]])
    for k, stop in enumerate(all_cuts[1:]):
        start = all_cuts[k]
        slc = seq[start - 1:stop - 1].seq.upper()
        for end in ends:
            if slc.startswith(end[0]) and slc.endswith(end[1]):
                if arguments.filter is not None:
                    if arguments.filter[0] <= len(slc) <= arguments.filter[1]:
                        fragments.append(len(slc))
                else:
                    fragments.append(len(slc))
                break
    if len(fragments) > 0:
        combo = "{0}-{1}".format(arguments.first, enzyme)
        print "\tCombo {}...".format(combo)
        cut_results[combo] = {
                'cuts': len(fragments),
                'cuts_per_bp': float(len(fragments)) / compare['seq_len'],
                'cut_length': numpy.mean(fragments),
                'cut_95_ci': 1.96 * numpy.std(fragments, ddof=1) / math.sqrt(len(fragments)),
                'min': min([i for i in fragments if i > 0]),
                'max': max(fragments)
            }
        return cut_results


def get_generic_end_results(arguments, enzyme, compare, all_cuts, cut_results):
    fragments = []
    for k, v in enumerate(all_cuts[1:]):
        slc = v - all_cuts[k]
        if arguments.filter is not None:
            if arguments.filter[0] <= slc <= arguments.filter[1]:
                fragments.append(slc)
        else:
            fragments.append(slc)
    if len(fragments) > 0:
        combo = "{0}-{1}".format(arguments.first, enzyme)
        print "\tCombo {}...".format(combo)
        cut_results[combo] = {
                'cuts': len(fragments),
                'cuts_per_bp': float(len(fragments)) / compare['seq_len'],
                'cut_length': numpy.mean(fragments),
                'cut_95_ci': 1.96 * numpy.std(fragments, ddof=1) / math.sqrt(len(fragments)),
                'min': min([i for i in fragments if i > 0]),
                'max': max(fragments)
            }
        return cut_results


def get_double_digest(arguments, single_cuts, seq):
    cut_results = {}
    base = single_cuts[arguments.first]
    print "\n\tDouble-cutting\n\t======================"
    for enzyme, compare in single_cuts.iteritems():
        # ensure we don't have any copy weirdness
        all_cuts = []
        all_cuts.extend(base['cuts'])
        if enzyme != arguments.first:
            all_cuts.extend(compare['cuts'])
            all_cuts = sorted(list(set(all_cuts)))
            if arguments.opposite_ends:
                cut_results = get_opposite_end_results(arguments, seq, base, enzyme, compare, all_cuts, cut_results)
            else:
                cut_results = get_generic_end_results(arguments, enzyme, compare, all_cuts, cut_results)
    return cut_results


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
        seq = SeqIO.read(open(fasta, 'rU'), 'fasta')
        single_cuts = get_single_digest(args, seq)
        results[fasta_name] = get_double_digest(args, single_cuts, seq)
        #pdb.set_trace()
    print "\nFile,Enzyme Combo,Cuts,Cuts per bp,Avg. fragment length,95CI cut length,Min(size),Max(size)"
    for fasta_name, values in results.iteritems():
        if values is not None:
            for combo_name, data in values.iteritems():
                print "{0},{1},{2[cuts]},{2[cuts_per_bp]},{2[cut_length]},{2[cut_95_ci]},{2[min]},{2[max]}".format(
                        fasta_name,
                        combo_name,
                        data)


if __name__ == '__main__':
    main()
