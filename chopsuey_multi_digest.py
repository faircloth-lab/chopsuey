#!/usr/bin/env python
# encoding: utf-8
"""
File: chopsuey.py
Author: Brant Faircloth

Created by Brant Faircloth on 06 August 2012 14:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Return summary data from a single virual digest of a directory of
targets using one or more restriction enzymes.

Usage:

    python chopsuey_single_digest.py /folder/with/some/fastas PstI SbfI BfaI NdeI AseI MseI

"""

import os
import sys
import argparse
import tempfile
import multiprocessing

from bx.seq import twobit

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
            description="""""")
    parser.add_argument(
            "--twobit",
            required=True,
            action=FullPaths,
            help="""A 2bit-formatted genome file"""
        )
    parser.add_argument(
            "--enzymes",
            required=True,
            type=str,
            nargs='+',
            help="""The enzyme(s) to chop with""",
        )
    parser.add_argument(
            "--output",
            required=True,
            action=FullPaths,
            help="""The output file to write results to.""",
        )
    parser.add_argument(
            "--min-frag-size",
            type=int,
            default=200,
            help="""The minimum frag size to accept.""",
        )
    parser.add_argument(
            "--max-frag-size",
            type=int,
            default=1000,
            help="""The maximum frag size to accept.""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of compute cores to use.""",
        )
    parser.add_argument(
            "--chunk-size",
            type=int,
            default=10000000,
            help="""Chunk size.""",
        )
    parser.add_argument(
        "--velvet",
        action="store_true",
        default=False,
        help="""Output data as velvet-style contigs""",
    )
    return parser.parse_args()


def chunk_scaffolds(genome, size):
    """Given a genome in many scaffolds, build temp files of `size` Mbp
    for easier querying"""
    print 'Chunking files into {0} bp...'.format(size)
    chromos = []
    tb = twobit.TwoBitFile(file(genome))
    # split target file into `options.size` (~10 Mbp) chunks
    fd, out = tempfile.mkstemp(suffix='.fasta')
    os.close(fd)
    temp = open(out, 'w')
    length = 0
    for seq in tb.keys():
        sequence = tb[seq][0:]
        length += len(sequence)
        # write it to the outfile
        temp.write('>{0}\n{1}\n'.format(seq, sequence))
        if length > size:
            # close current file and open a new one
            temp.close()
            # put tempfile name on stack
            chromos.append(out)
            # open a new temp file
            fd, out = tempfile.mkstemp(suffix='.fasta')
            os.close(fd)
            temp = open(out, 'w')
            # reset sequence length
            length = 0
    # ensure we close temp file
    if not temp.closed:
        temp.close()
    if out not in chromos:
        chromos.append(out)
    print "There are {} file chunks".format(len(chromos))
    return chromos


def worker(work):
    args, fasta, batch, idx = work
    # write new slices out to file
    fd, temp_name = tempfile.mkstemp(suffix='.fasta')
    # we need an cnt for all reads in chunk
    cnt = 0
    # there may be more than one sequence in a chunked file
    for seq in SeqIO.parse(open(fasta, 'rU'), 'fasta'):
        cuts = batch.search(seq.seq, linear=False)
        # merge all cuts across enzymes
        all_cuts = []
        for enzyme, sites in cuts.iteritems():
            all_cuts.extend(sites)
        # remove duplicates
        all = list(set(all_cuts))
        # sort
        all.sort()
        # add start and end of seq len to list
        all = [0] + all + [len(seq)]
        # create a list of tuples that are start and end coords
        slices = [(all[k], v) for k, v in enumerate(all[1:]) if args.min_frag_size <= v - all[k] < args.max_frag_size]
        for slice in slices:
            seqslice = seq[slice[0]:slice[1]]
            if args.velvet:
                seqslice.id = "NODE_{}{}|{}|{}-{}".format(idx, cnt, seqslice.id, slice[0], slice[1])
            else:
                seqslice.id = "{}|{}-{}".format(seqslice.id, slice[0], slice[1])
            seqslice.name = ""
            seqslice.description = ""
            os.write(fd, seqslice.format('fasta'))
            cnt += 1
    os.close(fd)
    sys.stdout.write(".")
    sys.stdout.flush()
    return temp_name


def main():
    args = get_args()
    print """
  ___  _  _   __  ____    ____  _  _  ____  _  _
 / __)/ )( \ /  \(  _ \  / ___)/ )( \(  __)( \/ )
( (__ ) __ ((  O )) __/  \___ \) \/ ( ) _)  )  /
 \___)\_)(_/ \__/(__)    (____/\____/(____)(__/

         Virtual Enzyme Single Digest
"""
    # chunk 2bit files
    chunks = chunk_scaffolds(args.twobit, args.chunk_size)
    # create a restriction batch
    batch = Restriction.RestrictionBatch(args.enzymes)
    work = [(args, chunk, batch, k) for k, chunk in enumerate(chunks)]
    if args.cores > 1:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores)
        results = pool.map(worker, work)
        # close the pool
        pool.close()
    else:
        results = map(worker, work)
    with open(args.output, 'w') as outfile:
        for fname in results:
            with open(fname) as infile:
                outfile.write(infile.read())

    # cleanup chunk files
    for fname in chunks:
        os.remove(fname)
    # cleanup temporary fasta results
    for fname in results:
        os.remove(fname)
    # done


if __name__ == '__main__':
    main()
