#!/usr/bin/env python

import argparse
import os
import numpy as np
from collections import Counter
from Bio import AlignIO
from argparse import RawDescriptionHelpFormatter


parser = argparse.ArgumentParser(description="Removes spurious characters from the beginning of a FASTA alignment.\n\n\
To trim multiple sequences, use the following command:\n\
\t\"for f in ./*.fa ; do ./pretrim -i $f ; done\"\n\n\
Or to trim in parallel, use:\n\
\t\"for f in ./A6.1_aligns/*.fa ; do echo ./pretrim.py -i $f ; done > run_pretrim.sh\"\n\
\t\"./parallel --jobs 8 < run_pretrim.sh\"",
                                 formatter_class=RawDescriptionHelpFormatter)

parser.add_argument("-align", "-i",
                    help="FASTA Alignment file",
                    required=True, type=str)

parser.add_argument("-threshold", "-t",
                    help="Minimum fraction of sites to consider as\
                          position as erroneous",
                    type=float, default=0.5)

parser.add_argument("-out", "-o",
                    help="Name of output file. Defaults to\
                    {input_alignment_name}.trim.fa",
                    required=False, type=str,
                    default=None)

args = parser.parse_args()

align = args.align
threshold = args.threshold
out = args.out

if threshold < 0 or threshold > 1:
    raise ValueError("threshold must be a fraction between 0 and 1")

del_empty = True


def read_align(align):
    alignment = AlignIO.read(align, "fasta")
    print("{} being processed".format(align))

    return alignment


def find_cutoff(alignment, threshold):
    to_cut = 0
    taxa = len(alignment)
    align_array = _align_to_array(alignment)
    empty_taxa = _count_empty(alignment)

    for col in align_array.T:
        sites = Counter(col)
        cnt_gap = sites[b"-"]
        gap_fract = float((cnt_gap-empty_taxa)/taxa)
        if gap_fract > threshold:
            to_cut += 1
        else:
            print("  Removing {} positions from beginning of {}".format(to_cut, align))
            trim_align = alignment[:, to_cut:]
            break
    return trim_align


def _align_to_array(alignment):
    align_array = np.array([list(rec) for rec in alignment],
                           np.character)
    return align_array


def _count_empty(alignment):
    empty_taxa = 0

    for record in alignment:
        sites = Counter(record.seq)
        seq_len = len(record)
        if sites["-"] == seq_len:
            print("\t{} consists of only gaps; ignoring".format(record.id))
            empty_taxa += 1
    return empty_taxa


def _del_empty(alignment):
    if del_empty:
        to_remove = []

        for record in alignment:
            sites = Counter(record.seq)
            seq_len = len(record)
            if sites["-"] == seq_len:
                print("\t{} consists of only gaps; deleting...".format(record.id))
                to_remove.append(record.id)

        if len(to_remove) > 0:
            print("  Sequences Removed")
            cleaned_records = [f for f in alignment if f.id not in to_remove]
            new_align = AlignIO.MultipleSeqAlignment(cleaned_records)

            return new_align

        else:
            print("  No sequences removed")
            return alignment

    else:
        return alignment



def write_out(trim_align, out):
    out_name = _get_outfile(out)
    AlignIO.write(trim_align, out_name, "fasta")


def _get_outfile(out):
    if out is None:
        out_name = "{}.trim.fa".format(os.path.splitext(align)[0])
    else:
        out_name = out
    return out_name


def main():
    alignment = read_align(align)
    trim_align = find_cutoff(alignment, threshold)

    trim_align2 = _del_empty(trim_align)
    
    write_out(trim_align2, out)


if __name__ == "__main__":
    main()
