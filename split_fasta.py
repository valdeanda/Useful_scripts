#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     split_fasta.py
# Purpose:  Split a fasta file according in almost equal parts
#           based on total base/residue count
#
# @uthor:      acph - dragopoot@gmail.com
# @mod:        vda  - valdeanda@utexas.edu
# Created:     Feb, 2018
# Copyright:   (c) acph 2018
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""Split a fasta file according in almost equal parts
based on total base/residue count.

Stores a numpy array that contains the lengths of the sequences in the file"""

import os
import argparse
import numpy as np
from subprocess import run, PIPE
from Bio import SeqIO


def arg_parse():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fastafile', type=str,
                        help='Fasta file to split')
    parser.add_argument('-p', '--parts', type=int, default=10,
                        help='Number of parts to slice the file [10]')
    args = parser.parse_args()
    return args


def count_entries(fname):
    """Count the number of entries in a fasta file.
    NOTE: Uses UNIX Grep.

    Keyword Arguments:
    fname -- Fasta file name.
    """
    cmd = ['grep', '-c', ">", fname]
    child = run(cmd, stdout=PIPE)
    count = int(child.stdout)
    return count


def create_array(fname, nentries, store=True):
    """Creates and stores the array with the length values for each sequence
    in the fasta file.

    Keyword Arguments:
    fname    -- Fasta file name
    nentries -- Number of entries in file
    """
    i = 0
    seqparser = SeqIO.parse(fname, 'fasta')
    lens = np.zeros(nentries, dtype=np.uint32)
    for seq in seqparser:
        lens[i] = len(seq)
        # if i % 1000000 == 0:
        #     print('{} ...'.format(i))
        i += 1
    # storing array
    if store:
        outfname = fname + '_lens.npy'
        np.save(outfname, lens)
    return lens


def split_fasta(fname, parts):
    """Split a fasta file into a number of parts base/residue wise.

    Keyword Arguments:
    fname -- Fasta file name
    parts -- Number of parts to split the file
    """
    seqparser = SeqIO.parse(fname, 'fasta')
    # load length array
    lenf = fname + '_lens.npy'
    if os.path.exists(lenf):
        print("Loading entries lengths...")
        lens = np.load(lenf)
        n = len(lens)
    else:
        print("Entries lengths array does not exist, creating:")
        print("   ...Calculing number of entries (using grep)")
        n = count_entries(fname)
        print("   ...Creating array!!!")
        lens = create_array(fname, n, store=True)
    # lenghts statistics
    lmean = lens.mean()
    lstd = lens.std()
    lsum = lens.sum()
    # limit of aminoacids per split (n=10)
    # the value is calculated substracting the average length and
    # 1 std length
    llimit = int((lsum-(lmean+lstd))/parts)
    # parts loop
    index = 0
    for i in range(parts):
        print('creating part{} of {} ...'.format(i+1, parts))
        outfname = fname + '.part{}-{}.fasta'.format(i+1, parts)
        outf = open(outfname, 'w')
        # initializing cummulative size
        cumsize = 0
        while cumsize < llimit:
            if index % 1000000 == 0:
                print('{} ...'.format(index))
            seq = next(seqparser, 'end')
            if type(seq) == str:
                break
            s = seq.seq
            d = seq.description
            t = '>{}\n{}\n'.format(d, s)
            outf.write(t)
            cumsize += len(s)
            index += 1
        # Close the file
        outf.close()


def main():
    """Main function

    """
    args = arg_parse()
    split_fasta(args.fastafile, args.parts)


if __name__ == '__main__':
    main()
