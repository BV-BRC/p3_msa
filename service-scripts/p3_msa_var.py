#!/usr/bin/env python

import argparse
import json
import sys


def run_msa_var(job_data, output_dir, tool_params):
    pass


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
