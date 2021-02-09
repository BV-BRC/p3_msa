#!/usr/bin/env python3
"""
Produces an MSA file consisting of columns of a given reference id
and a wiggle file of entropy for the reference id.
Removes gaps in the reference sequence.
Authors: Jacob S. Porter <jsporter@virginia.edu>
"""
import argparse
import datetime
import math
import os
import sys
from collections import Counter


def get_wiggle_from_MSA(
    MSA_file, seq_id, msa_output="msa.afa", wig_output="entropy.wig"
):
    """
    Produce an MSA file and a wiggle entropy file from a given MSA file
    and reference sequence

    Parameters
    ----------
    MSA_file: str
        The location of the MSA input file.
    seq_id: str
        The reference id for the sequence in the MSA file.
    msa_output: str
        The location to write the output MSA file to.
    wig_output: str
        The location to write the output WIG file to.

    Returns
    -------
    count: int
        A count of the length reference sequence after removing gaps.

    """
    count = 0
    order = []
    curr_id = ""
    seq = ""
    with open(MSA_file) as msa:
        for line in msa:
            if line.startswith(";") or line.startswith("#"):
                continue
            if line.startswith(">"):
                if curr_id:
                    order.append((curr_id, seq))
                curr_id = line[1:].strip()
                seq = ""
                continue
            seq += line.strip()
    if curr_id:
        order.append((curr_id, seq))
    stuff = [(i, rec) for i, rec in enumerate(order) if rec[0].startswith(seq_id)]
    if len(stuff) <= 0:
        return 0
    p, query = stuff[0]
    query_id, query_seq = query
    output_d = {q_id: "" for q_id, _ in order}
    del order[p]
    with open(msa_output, "w") as msa, open(wig_output, "w") as wig:
        print("variableStep chrom={}".format(query_id), file=wig)
        skip = 0
        for i, char in enumerate(query_seq):
            if char == "_" or char == "." or char == "~" or char == "-":
                skip += 1
                continue
            count += 1
            col = char
            output_d[query_id] += char
            for curr_id, curr_seq in order:
                col += curr_seq[i]
                output_d[curr_id] += curr_seq[i]
            if "n" in col or "N" in col:
                col = col.replace("n", "")
                col = col.replace("N", "")
            print("{} {}".format(i - skip, entropy(col)), file=wig)
        print(">{}".format(query_id), file=msa)
        print(output_d[query_id], file=msa)
        for curr_id, _ in order:
            print(">{}".format(curr_id), file=msa)
            print(output_d[curr_id], file=msa)
    return count


def entropy(seq, adj=100):
    """
    Calculate adjusted Shannon entropy from a sequence.

    Parameters
    ----------
    seq: str
        A sequence
    adj: int
        A constant to rescale Shannon entropy.

    Returns
    -------
    ent: int
        The adjusted entropy of the input.
        S = -100 * Sum (Pi * log2Pi)

    """
    counts = Counter(seq)
    ent = -adj * sum(
        [
            ((counts[char] / len(seq)) * math.log(counts[char] / len(seq), 2))
            for char in counts
        ]
    )
    return ent if ent > 0.0 else 0.0


def main():
    """Parse the arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser()
    parser.add_argument("MSA_file", help="Location of an MSA file.", type=str)
    parser.add_argument("seq_id", help="The sequence id of the reference.", type=str)
    parser.add_argument(
        "--msa_output",
        "-m",
        help="The location of the MSA output file.",
        type=str,
        default="msa.afa",
    )
    parser.add_argument(
        "--wig_output",
        "-w",
        help="The location of the WIG output file.",
        type=str,
        default="entropy.wig",
    )
    map_args = parser.parse_args()
    print("Started MSA_to_wiggle on {}".format(tick), file=sys.stderr)
    print(map_args, file=sys.stderr)
    sys.stderr.flush()
    count = get_wiggle_from_MSA(
        map_args.MSA_file, map_args.seq_id, map_args.msa_output, map_args.wig_output
    )
    print(
        "There were {} columns extracted.\nThe process took {} time.".format(
            count, datetime.datetime.now() - tick
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
