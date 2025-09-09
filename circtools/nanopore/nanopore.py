#! /usr/bin/env python3
#
# This file contains code that is derived from MIT-licensed software:
#   Original code by Morten T Venø <morten.veno@omiics.com>,
#   licensed under the MIT License.
#
# The MIT License:
# ---------------------------------------------------------------------
# See LICENSE_MIT for the full MIT license text.
# ---------------------------------------------------------------------
# Modifications and additions to this file are licensed under:
# The GNU General Public License (GPL), version 3 or later.
#
# Copyright (C) 2025 Tobias Jakobi
#
# See the LICENSE file in the root directory for the full terms of the GPL.
#
import re
from pathlib import Path
import subprocess

import csv
import gzip
import os
import sys

import pybedtools
import requests
from tqdm import tqdm

import yaml

import circ_module.circ_template
import circtools

def fix_path_for_docker(path: str):

    if is_running_in_docker() and os.path.isabs(path):
        path = str(os.path.join("/host_os/", path))
    return path


def is_running_in_docker():
    # Check for the presence of .dockerenv
    if os.path.exists('/.dockerenv'):
        return True

    return False

def is_writeable(directory):
    try:
        with open(os.path.join(directory, 'testfile'), 'w'):
            os.remove(os.path.join(directory, 'testfile'))
            pass
    except PermissionError:
        return False
    else:
        return True


def get_id_from_column_9(input_string: str, entity: str):
    splits = input_string.strip().split(';')

    feature_dict = {}

    for split in splits:
        item = split.strip().split(" ")
        if len(item) == 2:
            feature_dict[item[0]] = item[1].replace("\"", "")

    if entity in feature_dict:
        return feature_dict[entity]
    else:
        return None


def read_annotation_file(annotation_file, entity="exon"):
    """Reads a GTF file
    Will halt the program if file not accessible
    Returns a BedTool object only containing gene sections
    """
    """
    Reads a GTF file and outputs exons output used for the main script

    Download (needs to be the version__lift__version file as GTF, e.g. v37):
    https://www.gencodegenes.org/human/release_37lift37.html

    Expected output (each line one exon):
    chr1    11868   12227   DDX11L1 100     +       ENST00000456328.2_1     ENSG00000223972.5_4

    Args:
        annotation_file (str): Path to the GTF file.
    """
    try:
        file_handle = open(annotation_file)
    except PermissionError:
        message = ("Input file " + str(
            annotation_file) + " cannot be read, exiting.")
        sys.exit(message)
    else:

        with file_handle:
            line_iterator = iter(file_handle)
            bed_content = ""
            print("Start parsing GTF file")
            for line in line_iterator:
                # we skip any comment lines
                if line.startswith("#"):
                    continue

                # split up the annotation line
                columns = line.split('\t')

                if not (columns[2] == entity):
                    continue

                # we do not want any 0-length intervals -> bedtools segfault
                if int(columns[4]) - int(columns[3]) == 0:
                    continue

                entry = [
                    columns[0],
                    columns[3],
                    columns[4],
                    get_id_from_column_9(columns[8], "gene_name"),
                    "100",
                    columns[6],
                    get_id_from_column_9(columns[8], "transcript_id"),
                    get_id_from_column_9(columns[8], "gene_id")
                ]

                # concatenate lines to one string
                bed_content += '\t'.join(entry) + "\n"

        if not bed_content:
            exit(-1)

        # create a "virtual" BED file
        virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)

        return virtual_bed_file.sort()


def postprocess_ref_flat(refflat_csv: str):
    """
    Reads a refFlat file dumped by UCSC Genome Browser and outputs sorted
    exon output used for the main script

    refFlat schema:
    https://genome.ucsc.edu/cgi-bin/hgTables?hgta_doSchemaDb=hg38&hgta_doSchemaTable=refFlat

    Expected output:
    chr1    11868   12227   LOC102725121_exon_0_0_chr1_11869_f      0       +

    name = genename + exon # + 0 + chr + start+1 + strand (f or r)


    Args:
        file_path (str): Path to the input file.
        :param refflat_csv:
    """
    output_exons = ""
    output_genes = ""

    print("Creating refFlat-based exon files")
    try:
        with gzip.open(refflat_csv, mode='rt') as gz_file:
            csv_reader = csv.reader(gz_file, delimiter='\t')
            next(csv_reader, None)

            for row_number, row in enumerate(csv_reader, start=1):
                geneName, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds = row

                starts = exonStarts.split(',')
                stops = exonEnds.split(',')

                output_genes += "\t".join(
                    [chrom, txStart, txEnd, geneName, str(0),
                     strand]) + "\n"

                for exon_num in range(int(exonCount) - 1):
                    strand_tag = "f" if strand == "+" else "r"

                    name_tag = "_".join(
                        [geneName, "exon", str(exon_num), str(0), chrom,
                         str(int(starts[exon_num]) + 1), strand_tag])

                    output_genes += "\t".join(
                        [chrom, starts[exon_num], stops[exon_num], geneName,
                         str(0), strand]) + "\n"
                    output_exons += "\t".join(
                        [chrom, starts[exon_num], stops[exon_num], name_tag,
                         str(0), strand]) + "\n"

        virtual_bed_file = pybedtools.BedTool(output_exons,
                                              from_string=True)

        virtual_bed_file_genes = pybedtools.BedTool(output_genes,
                                                    from_string=True).sort()

        # the unique gene level file is different in regard to printing out all genes with comma instead of just choosing one
        # and omit the other co-optimal hits

        virtual_bed_file_genes = virtual_bed_file_genes.merge(s=True,
                                                              o=["distinct",
                                                                 "count_distinct",
                                                                 "distinct"],
                                                              c=[4, 4, 6])

        virtual_bed_file_sorted = virtual_bed_file.sort()

        virtual_bed_file_merged = virtual_bed_file_sorted.merge(s=True, o=[
            "collapse", "count", "distinct"], c=[4, 4, 6])

        file_base = refflat_csv.replace(".gz", "")

        with open(file_base + ".sort.bed", "w") as file:
            file.write(str(virtual_bed_file_sorted))

        with open(file_base + ".unique.bed", "w") as file:
            file.write(str(virtual_bed_file_genes))

        with open(file_base + ".merged.bed", "w") as file:
            file.write(str(virtual_bed_file_merged))

    except FileNotFoundError:
        print(f"Error: File '{refflat_csv}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


def postprocess_gencode(gencode_file: str):
    print("Creating GENCODE-based exon files")

    bed_file = read_annotation_file(gencode_file)

    file_base = gencode_file.replace(".gtf", "")

    with open(file_base + ".exon.bed", "w") as file:
        file.write(str(bed_file))

    virtual_bed_file_merged = bed_file.merge(s=True,
                                             o=["collapse",
                                                "count",
                                                "distinct"],
                                             c=[4, 4, 6])

    with open(file_base + ".exon.merge.bed", "w") as file:
        file.write(str(virtual_bed_file_merged))


def process_data(configuration: str, data_path: str):
    # build internal path from config name:

    with open(configuration, 'r') as config_file:
        config = (yaml.safe_load(config_file))

        full_data_path = os.path.join(data_path, config['dataset'])

        if not os.path.exists(data_path):
            os.mkdir(data_path)

        if is_writeable(data_path):
            print("Storing reference data in {}".format(data_path))

            # create folder, e.g. h19
            if not os.path.exists(full_data_path):
                os.makedirs(full_data_path)

            for item in config:

                if 'url' in config[item]:

                    url = config[item]['url']
                    file_name = os.path.join(full_data_path,
                                             config[item]['name'])

                    file_type = config[item]['type']

                    file_name_unzipped = file_name.replace("." + file_type,
                                                           "")

                    if not os.path.exists(file_name_unzipped):

                        with requests.get(url, stream=True) as r:
                            r.raise_for_status()
                            total_size = int(
                                r.headers.get('content-length', 0))

                            with open(file_name, 'wb') as f:
                                with tqdm(total=total_size, unit='B',
                                          unit_scale=True,
                                          desc="Downloading " +
                                               config[item][
                                                   'name']) as pbar:
                                    for chunk in r.iter_content(
                                            chunk_size=8192):
                                        f.write(chunk)
                                        pbar.update(len(chunk))

                        # most files need to be unpacked
                        if file_type == 'gz':
                            print("Unpacking.")
                            os.system("gzip -d " + file_name)
                            print("Done.")

                        # work on the gencode file
                        if 'postprocess' in config[item] and \
                                config[item]['postprocess'] == 'gencode':
                            postprocess_gencode(
                                gencode_file=file_name_unzipped)

                        # work with the refFlat file
                        if 'postprocess' in config[item] and \
                                config[item]['postprocess'] == 'refFlat':
                            postprocess_ref_flat(
                                refflat_csv=file_name_unzipped)

                    # file already exists
                    else:
                        print("Skipping, " + config[item]['name'] +
                              ", file already exists.")

class Nanopore(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version

        # mode
        self.run = self.cli_params.run
        self.download = self.cli_params.download
        self.check = self.cli_params.check

        # options
        self.threads =  self.cli_params.threads
        self.config = self.cli_params.config
        self.sample_name = self.cli_params.sample
        self.sample_path = self.cli_params.sample
        self.reference_path = self.cli_params.reference_path
        # fix for docker if default argument is passed
        self.reference_path = fix_path_for_docker(self.reference_path)
        self.keep_temp = self.cli_params.keep_temp
        self.threads = self.cli_params.threads
        self.dry_run = self.cli_params.dry_run
        self.output_path = self.cli_params.output_path

        # set path to library directory
        self.script_path = circtools.__path__[0] + "/nanopore/"
        self.config_location = circtools.__path__[0] + "/nanopore/config"


    def module_name(self):
        """Return a string representing the name of the module."""
        return self.program_name

    def run_module(self):

        if self.cli_params.run:
            self.run_analysis()
        elif self.cli_params.download:
            self.run_download()
        elif self.cli_params.check:
            self.run_check(verbose=True)

    def run_download(self):

        if not self.config:
            print("No configuration file specified.")
            exit(-1)

        final_path = os.path.join(self.config_location, self.config+".yml")

        if not os.path.isfile(final_path):
            print("Configuration file {} not accessible.".format(final_path))
            exit(-1)
        else:
            process_data(configuration=str(final_path),
                         data_path=self.reference_path)


    def run_check(self, verbose = False):

        tools = {"bedtools": "bedtools",
                 "NanoFilt": "NanoFilt -h",
                # "pblat": "pblat",
                 "samtools": "samtools --help"}

        failed = False
        warn = False

        for tool in tools:
            if verbose:
                print("Checking for {}".format(tool))
            try:
                output = subprocess.run(tools[tool].split(" "),
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, timeout=5)
                if output.returncode != 0:
                    print("\tThere might be a problem with {}".format(tool))
                    warn = True
                    failed = True
            except FileNotFoundError:
                print("\tUnable to find {}!".format(tool))
                failed = True
            except Exception as e:
                print(e)

        if warn:
            print()
            print("WARNING: There might be a problem "
                  "with some of the software installed")
            print()

        if failed:
            print()
            print("ERROR: Some of the required software is missing!")
            print()
            sys.exit(1)
        else:
            if verbose:
                print()
                print("All of the expected software requirements are present!")
                print()


    def run_analysis(self):

        # always do pre-flight check before analysis
        self.run_check(verbose=False)

        if not self.sample_path:
            print("No sample path specified.")
            exit(-1)

        # Check if sample_path exists
        if not os.path.exists(self.sample_path):
            print(
                "Error: '{}' does not exists! Please make sure that the path is written correctly".format(
                    self.sample_path))
            exit(-1)
        else:
            self.sample_path = os.path.dirname(os.path.abspath(self.cli_params.sample))

        if not self.sample_name:
            print("No sample file specified.")
            exit(-1)

        if not self.config:
            print("No genome build configuration specified.")
            exit(-1)

        if not self.output_path:
            print("No output path specified.")
            exit(-1)

        # Check if sample exists
        if not os.path.exists(self.sample_name):
            print("Error: Sample file '{}' does not exist!".format(self.sample_name))
            exit(-1)

        # Check if sample ends with ".fq.gz"
        # if self.sample_name and not self.sample_name.endswith(".fq.gz"):
        #     print("Error: Sample file '{}' does not end with '.fq.gz'".format(
        #         self.sample_name))
        #     exit(-1)

        # Check if reference_path exists
        if not os.path.exists(self.reference_path):
            print(
                "Error: '{}' does not exists! Please make sure that the path is written correctly".format(
                    self.reference_path))
            exit(-1)
        else:
            self.reference_path = os.path.abspath(self.reference_path)

        # Check if reference_path is a directory:
        if not os.path.isdir(self.reference_path):
            print(
                "Error: '{}' is not a directory! Please provide a directory that contains the reference data.".format(
                    self.reference_path))
            exit(-1)

        # Check if reference_path contains genome in the end
        if os.path.basename(self.reference_path) == self.config:
            print(
                "Error: Looks like you have provided the direct path to the genome directory '{}'. Please simply provide data directory '{}'.".format(
                    self.reference_path, os.path.dirname(self.reference_path)
                ))
            exit(-1)

        # Check if reference_path + genome exists
        if self.reference_path and self.config and not os.path.exists(os.path.join(self.reference_path, self.config)):
            print(
                "Error: The genome build {} in the directory {} does not exist!".format(
                    self.config, self.reference_path))
            print(
                "Note: You might need to download genome data first using the download command.")
            exit(-1)

        # Prepare sample_path and sample_name

        # sample_name = re.sub(r'\.\w+\.\w+$', "",self.sample_name)

        #ext = re.sub(r'\.\w+\.\w+$', "",self.sample_name)

        # Create a regex pattern with groups for name and age
        pattern = r'([a-zA-Z0-9_-]+)\.(\w+\.\w+)$'

        # Find a match in the text
        match = re.search(pattern, self.sample_name)

        if match:
            # Access captured groups
            sample_ext = match.group(2)
            sample_name = match.group(1)


        # sample_name = os.path.basename(self.sample_name).replace('.fq.gz', '')
        script_path = str(Path(os.path.expanduser(self.script_path)).resolve())
        output_path = str(Path(self.output_path).resolve())

        final_path = os.path.join(self.config_location, self.config+".yml")


        # Check if the script-path exists
        if not os.path.exists(script_path):
            print(
                "Error: The script path '{}' does not exists! Please point to the location of the scripts directory in the long_read_circRNA path using --script-path".format(
                    self.script_path))
            exit(-1)

        # Check if the required scripts are in the script script-path
        target_scripts = ["blat_nanopore_v6.0.sh",
                          "novel_exons_and_alternative_usage_v8.0.sh"]
        for target_script in target_scripts:
            if not os.path.exists(os.path.join(self.script_path, target_script)):
                print(
                    "Error: '{}' script is not found in the script path '{}'. Are you sure the provided script-path is correct?".format(
                        target_script, self.script_path))
                exit(-1)

        # Check if the output path exists
        if not os.path.exists(output_path):
            os.mkdir(output_path)
            print("Output_path '{}' does not exist, "
                  "creating directory.".format(self.output_path))

        #output_path = os.path.join(self.output_path, self.sample_name)

        print("Starting process with the following settings")
        print("Sample name: {}".format(sample_name))
        print("Reference path: {}".format(self.reference_path))
        print("Sample path: {}".format(self.sample_path))
        print("Genome build: {}".format(self.config))
        print("Script-path: {}".format(self.script_path))
        print("Output-path: {}".format(self.output_path))
        print()

        original_directory = os.getcwd()

        if not self.dry_run:
            # Main process for circRNA detection
            print( os.path.join(self.script_path, "blat_nanopore_v6.0.sh"))
            subprocess.run(
                ["bash", os.path.join(self.script_path, "blat_nanopore_v6.0.sh"),
                 self.sample_path,
                 sample_name,
                 self.config,
                 self.reference_path,
                 self.script_path,
                 self.output_path,
                 str(self.threads),
                 sample_ext],
                 )
            print("")

            print("circRNA detection has finished")

            with open(final_path, 'r') as config_file:
                build_config = (yaml.safe_load(config_file))

                if len(build_config) == 7:

                    print("Starting the novel exon and alternative usage script")

                    os.chdir(original_directory)

                    keep_check = {True: "yes", False: "no"}

                    keep_temp = keep_check[self.keep_temp]

                    subprocess.run(["bash",
                                    os.path.join(self.script_path,
                                    "novel_exons_and_alternative_usage_v8.0.sh"),
                                    sample_name,
                                    self.config,
                                    self.reference_path,
                                    self.script_path,
                                    self.output_path,
                                    keep_temp,
                                    str(self.threads)
                                    ]
                                   )
                    print("")

                else:
                    print("Skipping novel exon and alternative usage script.")
                    print("Only supported for human and mouse genome.")

            print("Long_read_circRNA has finished!")
        else:
            print("Dry run is complete")
