#!/usr/bin/env python3
"""
Produces an MSA file consisting of columns of a given reference id
and a wiggle file of entropy for the reference id.
Removes gaps in the reference sequence.
TODO: this could be done in parallel.
TODO: this could be done with seek on a binary file to reduce memory requirements.
Authors: Jacob S. Porter <jsporter@virginia.edu>
"""
import argparse
import datetime
import math
import os
import sys
from collections import Counter
from math import e, log

import numpy as np
from tqdm import tqdm


def get_wiggle_from_MSA(
    MSA_file, seq_id, wig_output="entropy.wig", temp_output=".msa.temp"
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
    # output_d = {q_id: "" for q_id, _ in order}
    del order[p]
    with open(temp_output, "w") as msa, open(wig_output, "w") as wig:
        print("variableStep chrom={}".format(query_id), file=wig)
        skip = 0
        print(
            "Length of reference: {}".format(len(query_seq.replace("-", ""))),
            file=sys.stderr,
        )
        sys.stderr.flush()
        print(">{}".format(query_id), file=msa)
        for curr_id, _ in order:
            print(">{}".format(curr_id), file=msa)
        # output_reference = []
        for i, char in enumerate(tqdm(query_seq)):
            if i % 5000 == 1:
                wig.flush()
                msa.flush()
            if char == "_" or char == "." or char == "~" or char == "-":
                skip += 1
                continue
            count += 1
            col = [char]
            # output_d[query_id].append(char)
            msa.write(char)
            for curr_id, curr_seq in order:
                if curr_seq[i].upper() != "N":
                    col.append(curr_seq[i])
                msa.write(curr_seq[i])
                # output_d[curr_id] += curr_seq[i]
            print("{} {}".format(i - skip, entropy2(col)), file=wig)
        # print(">{}".format(query_id), file=msa)
        # print("".join(output_d[query_id]), file=msa)
        # for curr_id, _ in order:
        # print(">{}".format(curr_id), file=msa)
        #     print("".join(output_d[curr_id]), file=msa)
    return count


def convertFasta(msa_output="msa.afa", temp_output=".msa.temp"):
    with open(temp_output) as temp, open(msa_output, "w") as msa:
        ids = []
        seq = None
        for line in temp:
            if line.startswith(">"):
                ids.append(line.strip())
            else:
                seq = line.strip()
        for i, record_id in tqdm(enumerate(ids)):
            print(record_id, file=msa)
            for j, char in enumerate(seq):
                if j % len(ids) == i:
                    msa.write(char)
            msa.write("\n")
    if os.path.exists(temp_output):
        os.remove(temp_output)


def entropy1(seq, adj=100):
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


def entropy2(labels, base=None, adj=100):
    """
    Computes entropy of label distribution.
    Modified from:
    https://stackoverflow.com/questions/15450192/fastest-way-to-compute-entropy-in-python
    """
    labels = list(labels)
    n_labels = len(labels)
    if n_labels <= 1:
        return 0.0
    _, counts = np.unique(labels, return_counts=True)
    probs = counts / n_labels
    n_classes = np.count_nonzero(probs)
    if n_classes <= 1:
        return 0.0
    ent = 0.0
    # Compute entropy
    base = e if base is None else base
    for i in probs:
        ent -= i * log(i, base)
    return adj * ent


def main():
    """Parse the arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    parser_wiggle = subparsers.add_parser("wiggle", help="Produce wiggle and afa file.")
    parser_wiggle.add_argument("MSA_file", help="Location of an MSA file.", type=str)
    parser_wiggle.add_argument(
        "seq_id", help="The sequence id of the reference.", type=str
    )
    parser_wiggle.add_argument(
        "--msa_output",
        "-m",
        help="The location of the MSA output file.",
        type=str,
        default="msa.afa",
    )
    parser_wiggle.add_argument(
        "--wig_output",
        "-w",
        help="The location of the WIG output file.",
        type=str,
        default="entropy.wig",
    )
    parser_fasta = subparsers.add_parser(
        "fasta", help="Produce afa file from intermediate temp file."
    )
    parser_fasta.add_argument(
        "temp_file", help="The location for the intermediate temp file.", type=str
    )
    parser_fasta.add_argument(
        "--msa_output",
        "-m",
        help="The location of the MSA output file.",
        type=str,
        default="msa.afa",
    )
    map_args = parser.parse_args()
    mode = map_args.mode
    if mode == "wiggle":
        print("Started MSA_to_wiggle on {}".format(tick), file=sys.stderr)
        print(map_args, file=sys.stderr)
        sys.stderr.flush()
        count = get_wiggle_from_MSA(
            map_args.MSA_file, map_args.seq_id, map_args.wig_output
        )
        convertFasta(map_args.msa_output)
        print("There were {} columns extracted.".format(count), file=sys.stderr)
    elif mode == "fasta":
        convertFasta(map_args.msa_output, map_args.temp_file)
    else:
        print("No such mode: {}".format(mode), file=sys.stderr)
    print(
        "The process took {} time.".format(datetime.datetime.now() - tick),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
