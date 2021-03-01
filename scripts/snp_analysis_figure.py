#!/usr/bin/env python
"""
Create figures for variance analysis.
Jacob S. Porter
"""
import argparse
import csv

import matplotlib.pyplot as plt

plt.style.use("ggplot")


def get_figures(table_location, output_prefix):
    position = []
    score = []
    counts = []
    with open(table_location, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            try:
                position.append(int(row["Position"]))
                score.append(int(row["Score"]))
                counts.append(int(row["NumberOfSequence"]))
            except ValueError:
                continue
    fig, ax1 = plt.subplots()
    color = "tab:blue"
    ax1.set_xlabel("Position")
    ax1.set_ylabel("Score", color=color)
    ax1.bar(position, score, color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    ax2 = ax1.twinx()
    color = "tab:red"
    ax2.set_ylabel("Num of Seqs.", color=color)
    ax2.plot(position, counts, color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    plt.title("Position versus Score and Number of Sequences")
    fig.tight_layout()
    # plt.bar(position, score, color="blue")
    # plt.xlabel("Position")
    # plt.ylabel("Score")
    plt.savefig(output_prefix + ".svg", format="svg")
    plt.savefig(output_prefix + ".png", format="png")
    return len(position)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("table", help="Location of SNP entropy table.", type=str)
    parser.add_argument(
        "output_prefix",
        help="The location and prefix for saving the figures.",
        type=str,
    )
    map_args = parser.parse_args()
    get_figures(map_args.table, map_args.output_prefix)


if __name__ == "__main__":
    main()
