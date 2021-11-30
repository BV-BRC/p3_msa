#!/usr/bin/env python3
"""
Convert aligned fasta file into other formats.
Jacob S. Porter
"""
import argparse
import datetime
import sys
from multiprocessing import Process

from Bio import AlignIO

INFORMAT = "fasta"
OUTFORMATS = [
    ("clustal", "", "aln"),
    ("nexus", "", "nexus"),
    ("phylip-relaxed", "", "phy"),
    ("pir", "", "pir"),
    # ("phylip", "", "phylip"), # phylip truncates seq ids to 10 characters and needs them to be unique.
    # ("phylip-sequential", ".seq", "phylip"),
]


def convert_file(in_file, out_file_prefix, in_format=INFORMAT, molecule_type=None):
    print(in_file, out_file_prefix, file=sys.stderr)
    p_list = []
    count = 0
    for form, form_txt, ending in OUTFORMATS:
        print("{} {}".format(form, ending), file=sys.stderr)
        p = Process(
            target=AlignIO.convert,
            args=(
                in_file,
                in_format,
                "{}{}.{}".format(out_file_prefix, form_txt, ending),
                form,
                molecule_type,
            ),
        )
        p.start()
        p_list.append(p)
    for p in p_list:
        p.join()
        count += 1
    return count


def main():
    now = datetime.datetime.now()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("in_file", help="Location of input file.", type=str)
    parser.add_argument("out_file_prefix", help="Prefix of output file.", type=str)
    parser.add_argument(
        "molecule_type",
        help="Molecule type to apply: “DNA”, “RNA” or “protein”.",
        type=str,
        choices=["DNA", "protein", "RNA"],
    )
    parser.add_argument(
        "--in_format", help="The input file format.", type=str, default=INFORMAT
    )
    map_args = parser.parse_args()
    print("Converting alignment file to other formats.", file=sys.stderr)
    print(map_args, file=sys.stderr)
    count = convert_file(
        map_args.in_file,
        map_args.out_file_prefix,
        map_args.in_format,
        map_args.molecule_type,
    )
    print(
        "Converted {} files in time: {}".format(count, datetime.datetime.now() - now),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
