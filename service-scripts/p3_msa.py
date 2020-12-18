#!/usr/bin/env python3
"""
Run the MSA service and analyze variation.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import requests

import fqutil_api as patric_api

# until / unless we get around to pip'ing. enable local and installed versions to work
# if os.path.exists("../lib/fqutil_api.py"):
#    sys.path.insert(1, '../lib')

# The max number of characters to check when checking a fasta file for AA or NT
LINE_LEN = 1000


def get_feature_sequences(feature_grp_path):
    target_file = "feature_stuff.txt"
    feature_grp_path = "/jsporter@patricbrc.org/home/MSA/feature_group_eg/AlignMe1"
    feature_url = "https://p3.theseed.org/services/data_api/genome_feature/?in(feature_id,FeatureGroup({}))&limit(25000)".format(
        feature_grp_path)
    headers = {"accept": "application/json"}
    req = requests.Request('GET', feature_url, headers=headers)
    print(feature_url)
    patric_api.authenticateByEnv(req)
    prepared = req.prepare()
    my_sess = requests.Session()
    response = my_sess.send(prepared)
    handle = open(target_file, 'wb')
    if not response.ok:
        sys.stderr.write("API not responding. Please try again later.\n")
        sys.exit(2)
    else:
        for block in response.iter_content(1024):
            handle.write(block)


def check_nt(file_path):
    """Check if the file_path is a NT fasta file.  Else, it is AA."""
    with open(file_path, "r") as fasta_file:
        test_line = ""
        for line in fasta_file:
            if not line.startswith(">") and not line.startswith(
                    "#") and not line.startswith(";"):
                test_line += line.strip().upper()[:LINE_LEN]
            if len(test_line) >= LINE_LEN:
                break
        test_set = set(test_line)
        cmp_set = set("ACTGN")
        if len(test_set) <= 5 and len(
                test_set.intersection(cmp_set)) == len(test_set):
            return True
        return False


def run_msa(job_data, output_dir, tool_params):
    basenames = set()
    count = 1
    for file_object in job_data["fasta_files"]:
        my_output_dir = os.path.join(output_dir, str(count))
        os.mkdir(my_output_dir)
        os.symlink(file_object["file"],
                   os.path.join(my_output_dir, "input.fasta"))
        var_cmd = [
            "/homes/jsporter/p3_msa/p3_msa/service-scripts/web_flu_snp_analysis.pl",
            "-r", my_output_dir
        ]
        nucl = check_nt(file_object["file"])
        if nucl:
            var_cmd += ["-n"]
        basename = os.path.basename(file_object["file"])
        base_count = 2
        basename_temp = basename
        while basename_temp in basenames:
            basename_temp = basename + "_" + str(base_count)
            base_count += 1
        basenames.add(basename_temp)
        basename = basename_temp.split(".")
        basename = ".".join(basename[0:len(basename) - 1])
        subprocess.check_call(var_cmd)
        # os.unlink(os.path.join(output_dir, "input.fasta"))
        # Save: output.aln, output.afa, cons.fasta, foma.table
        # Ignore: gaf.log, main.log, output.zip, snp.out
        shutil.move(os.path.join(my_output_dir, "output.aln"),
                    os.path.join(output_dir, "{}.aln".format(basename)))
        shutil.move(os.path.join(my_output_dir, "output.afa"),
                    os.path.join(output_dir, "{}.afa".format(basename)))
        shutil.move(os.path.join(my_output_dir, "cons.fasta"),
                    os.path.join(output_dir, "{}.cons.fasta".format(basename)))
        shutil.move(os.path.join(my_output_dir, "foma.table"),
                    os.path.join(output_dir, "{}.foma.tsv".format(basename)))
        count += 1


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
    get_feature_sequences("")
    # run_msa(job_data, map_args.o, tool_params)


if __name__ == "__main__":
    main()
