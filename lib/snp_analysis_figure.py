#!/usr/bin/env python
import argparse
import csv

import matplotlib.pyplot as plt
"""
Create figures for variance analysis.
Jacob S. Porter
"""
plt.style.use('ggplot')


def get_figures(table_location, output_prefix):
    position = []
    score = []
    with open(table_location, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            try:
                position.append(int(row["Position"]))
                score.append(int(row["Score"]))
            except ValueError:
                continue
    plt.bar(position, score, color='blue')
    plt.xlabel("Position")
    plt.ylabel("Score")
    plt.title("Position versus Score")
    plt.savefig(output_prefix + ".svg", format="svg")
    plt.savefig(output_prefix + ".png", format="png")
    return len(position)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('table',
                        help='Location of SNP entropy table.',
                        type=str)
    parser.add_argument('output_prefix',
                        help='The location and prefix for saving the figures.',
                        type=str)
    map_args = parser.parse_args()
    get_figures(map_args.table, map_args.output_prefix)


if __name__ == "__main__":
    main()
