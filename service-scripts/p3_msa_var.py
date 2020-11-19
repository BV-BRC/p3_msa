#!/usr/bin/env python3

import argparse
import json
import os
import shutil
import subprocess
import sys

# The max number of characters to check when checking a fasta file for AA or NT
LINE_LEN = 1000


def check_nt(file_path):
    """Check if the file_path is a NT fasta file.  Else, it is AA."""
    with open(file_path, "r") as fasta_file:
        test_line = ""
        for line in fasta_file:
            if not line.startswith(">"):
                test_line += line.strip().upper()[:LINE_LEN]
            if len(test_line) >= LINE_LEN:
                break
        test_set = set(test_line)
        cmp_set = set("ACTGN")
        if len(test_set) <= 5 and len(
                test_set.intersection(cmp_set)) == len(test_set):
            return True
        return False


def run_msa_var(job_data, output_dir, tool_params):
    for file_object in job_data["fasta_files"]:
        os.symlink(file_object["file"], os.path.join(output_dir,
                                                     "input.fasta"))
        var_cmd = [
            "/homes/jsporter/p3_msa_var/p3_msa_var/service-scripts/web_flu_snp_analysis.pl",
            "-r", output_dir
        ]
        if check_nt(file_object["file"]):
            var_cmd += ["-n"]
        subprocess.check_call(var_cmd)
        os.unlink(os.path.join(output_dir, "input.fasta"))
        # output.aln, output.afa, cons.fasta, foma.table
        shutil.move(os.path.join(output_dir, "foma.table"),
                    os.path.join(output_dir, "foma.tsv"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--jfile', help='json file for job', required=True)
    parser.add_argument(
        '--sstring',
        help='json server string specifying api {"data_api":"url"}',
        required=True,
        default=None)
    parser.add_argument(
        '-p',
        help='JSON formatted parameter list for info keyed to program',
        default='{}',
        required=False)
    parser.add_argument(
        '-o',
        help='output directory. Defaults to current directory.',
        required=False,
        default="./")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    map_args = parser.parse_args()
    with open(map_args.jfile, 'r') as job_handle:
        job_data = json.load(job_handle)
    server_info = json.loads(map_args.sstring)
    for k, d in server_info.items():
        job_data[k] = d
    job_data["output_path"] = map_args.o
    try:
        tool_params = json.loads(map_args.p)
    except json.decoder.JSONDecodeError:
        tool_params = {}
    print("Parameters: {}".format(tool_params), file=sys.stdout)
    run_msa_var(job_data, map_args.o, tool_params)


if __name__ == "__main__":
    main()
