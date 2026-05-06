#!/usr/bin/env python3

import argparse
import os
import re
import subprocess

r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "circtools_exon_wrapper.R")


def main():
    parser = argparse.ArgumentParser(description="Exon usage analysis (ballgown + edgeR + CircTest)")
    parser.add_argument("-d", "--detect-dir",        required=True,  dest="detect_dir")
    parser.add_argument("-r", "--replicates",        required=True,  dest="replicates")
    parser.add_argument("-c", "--condition-list",    required=True,  dest="condition_list")
    parser.add_argument("-l", "--condition-columns", required=True,  dest="condition_columns")
    parser.add_argument("-g", "--grouping",          required=True,  dest="grouping")
    parser.add_argument("-o", "--output-directory",  required=True,  dest="output_directory")
    parser.add_argument("-p", "--output-prefix",     required=True,  dest="output_prefix")
    parser.add_argument("-b", "--ballgown-data",     required=True,  dest="ballgown_data")
    parser.add_argument("-G", "--gtf-file",          required=True,  dest="gtf_file")
    parser.add_argument("-T", "--circtest-file",     required=True,  dest="circtest_file")
    parser.add_argument("--has-header",              action="store_true", dest="has_header")
    parser.add_argument("-s", "--species",           required=True,  dest="species",
                        choices=["hs", "rn", "mm", "ss"])
    cli_params = parser.parse_args()

    if not os.access(cli_params.output_directory, os.W_OK):
        print("Output directory %s not writable." % cli_params.output_directory)
        exit(-1)

    if not os.path.exists(cli_params.output_directory):
        os.makedirs(cli_params.output_directory)

    if not os.path.exists(cli_params.detect_dir):
        print("DCC/detect data directory %s does not exist (or is not accessible)." % cli_params.detect_dir)
        exit(-1)

    for f in ["CircRNACount", "LinearCount", "CircCoordinates"]:
        path = cli_params.detect_dir + f
        if not os.path.exists(path):
            print("Input file %s not found." % path)
            exit(-1)

    if cli_params.ballgown_data and not os.path.exists(cli_params.ballgown_data):
        print("Ballgown data directory %s does not exist (or is not accessible)." % cli_params.ballgown_data)
        exit(-1)

    for f in [cli_params.gtf_file, cli_params.circtest_file]:
        if not os.path.exists(f):
            print("Input file %s not found." % f)
            exit(-1)

    if len(cli_params.condition_list.split(",")) < 2:
        print("Error: Length of parameter list specified via -c is < 2.")
    if len(cli_params.condition_columns.split(",")) < 2:
        print("Error: Length of parameter list specified via -l is < 2.")
    if len(cli_params.grouping.split(",")) < 2:
        print("Error: Length of parameter list specified via -g is < 2.")

    for column in cli_params.condition_columns.split(","):
        try:
            int(column)
        except ValueError:
            print("Error: column %s is no valid column index." % str(column))
            exit(-1)

    for column in cli_params.grouping.split(","):
        try:
            int(column)
        except ValueError:
            print("Error: group %s is no valid group index." % str(column))
            exit(-1)

    command = "Rscript"
    r_location = subprocess.check_output(['which', command], universal_newlines=True,
                                         stderr=subprocess.STDOUT).split('\n')[0]
    r_version = subprocess.check_output([command, '--version'], universal_newlines=True,
                                        stderr=subprocess.STDOUT)
    m = re.search(r'(\d+\.\d+\.\d+)', r_version)
    r_version = m.group(0)
    print("Using R version %s [%s]" % (r_version, r_location))

    args = [
        cli_params.detect_dir,
        cli_params.replicates,
        cli_params.condition_list,
        cli_params.condition_columns,
        cli_params.grouping,
        cli_params.output_directory + "/" + cli_params.output_prefix,
        cli_params.ballgown_data,
        cli_params.gtf_file,
        cli_params.circtest_file,
        cli_params.has_header,
        cli_params.species,
    ]

    cmd = ["Rscript", r_script] + [str(e) for e in args]
    subprocess.call(cmd)


if __name__ == "__main__":
    main()