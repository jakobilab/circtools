#! /usr/bin/env python3

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import circ_module.circ_template

import os
import sys
import string
import random
import itertools
import subprocess

import pybedtools
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram

class Conservation(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version
        self.temp_dir = self.cli_params.global_temp_dir
        self.gtf_file = self.cli_params.gtf_file
        self.fasta_file = self.cli_params.fasta_file
        self.detect_dir = self.cli_params.detect_dir
        self.output_dir = self.cli_params.output_dir
        self.organism = self.cli_params.organism
        self.id_list = self.cli_params.id_list
        self.experiment_title = self.cli_params.experiment_title
        self.input_circRNA = self.cli_params.sequence_file
        
        # gene_list file argument
        if (self.cli_params.gene_list):
            self.gene_list = self.cli_params.gene_list
        elif (self.cli_params.gene_list_file):
            gene_file = [line.rstrip() for line in open(self.cli_params.gene_list_file[0])]
            self.gene_list = gene_file
        else:
            print("Need to provide gene list by either -G or -GL options!")
            exit(-1)

        if self.id_list and self.gene_list:
            print("Please specify either host genes via -G/-GL or circRNA IDs via -i.")
            sys.exit(-1)

        self.other_blast_db = "nt"

        if self.organism not in ["mm", "hs", "ss", "rn", "cl"]:
            print("Please provide valid species. Options available: Mouse, Human, Pig, Rat, Dog")
            sys.exit(-1)

    def module_name(self):
        """Return a string representing the name of the module."""
        return self.program_name

    # Register an handler for the timeout
    def handler(self, signum, frame):
        raise Exception("Maximum execution time for remote BLAST reached. Please try again later.")

    @staticmethod
    def read_annotation_file(annotation_file, entity="gene", string=False):
        """Reads a GTF file
        Will halt the program if file not accessible
        Returns a BedTool object only containing gene sections
        """

        try:
            file_handle = open(annotation_file)
        except PermissionError:
            message = ("Input file " + str(annotation_file) + " cannot be read, exiting.")
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

                    # extract chromosome, start, stop, score(0), name and strand
                    # we hopefully have a gene name now and use this one for the entry
                    
                    # added by Shubhada (to fetch the gene names)
                    s = str(columns[8])
                    news = s.strip("\n")[:-1].replace("; ", ";")          # removing trailing ; to form dictionary in next step
                    temp = [x for x in news.replace("\"", "").split(";")]
                    temp_keys = [x.split(" ")[0] for x in temp]
                    temp_values = ["_".join(x.split(" ")[1:]) for x in temp]
                    gene_dict = dict(zip(temp_keys,temp_values))
                    if ("gene_name" in gene_dict.keys()):
                        gene_name = gene_dict["gene_name"]
                    else:
                        gene_name = "name"
                    entry = [
                        columns[0],
                        columns[3],
                        columns[4],
                        gene_name,
                        str(0),
                        columns[6]
                        #columns[1]              # flag ensemble/havana
                    ]

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"

            if not bed_content:
                exit(-1)

            if string:
                return bed_content
            else:
                return bed_content

    def run_module(self):

        if self.id_list and os.access(self.id_list[0], os.R_OK):
            print("Detected supplied circRNA ID file.")
            with open(self.id_list[0]) as f:
                lines = f.read().splitlines()
            self.id_list = lines

        # let's first check if the temporary directory exists
        if not (os.access(self.temp_dir, os.W_OK)):
            print("Temporary directory %s not writable." % self.temp_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        # let's first check if the temporary directory exists
        if not (os.access(self.output_dir, os.W_OK)):
            print("Output directory %s not writable." % self.output_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        circ_rna_number = 0

        # call the read_annotation_file and store exons in both bed and bedtools format for linear and circRNA
        exons_bed = self.read_annotation_file(self.gtf_file, entity="exon")
        exons_bed_list = [x.split("\t") for x in exons_bed.strip().split("\n")]
        # create a "virtual" BED file for circular RNA bedtools intersect
        virtual_bed_file = pybedtools.BedTool(exons_bed, from_string=True)
        print("Start merging GTF file outside the function")
        # we trust that bedtools >= 2.27 is installed. Otherwise this merge will probably fail
        exons = virtual_bed_file.sort().merge(s=True,  # strand specific
                                                 c="4,5,6",  # copy columns 5 & 6
                                                 o="distinct,distinct,distinct")  # group
        #print(exons_bed_list[:5])
        exons_bed_list = [x.split("\t") for x in str(exons).splitlines()]

        flanking_exon_cache = {}
        all_exons_circle = {} 
        if self.detect_dir:
            with open(self.detect_dir) as fp:
                
                for line in fp:

                    # make sure we remove the header
                    if line.startswith('Chr\t'):
                        continue

                    line = line.rstrip()
                    current_line = line.split('\t')
                    if current_line[3] == "not_annotated":
                        continue

                    if self.gene_list and not self.id_list and current_line[3] not in self.gene_list:
                        continue
                        
                    sep = "_"
                    name = sep.join([current_line[3],
                                        current_line[0],
                                        current_line[1],
                                        current_line[2],
                                        current_line[5]])

                    if self.id_list and not self.gene_list and name not in self.id_list:
                        continue

                    circrna_length = int(current_line[2]) - int(current_line[1])

                    sep = "\t"
                    bed_string = sep.join([current_line[0],
                                            current_line[1],
                                            current_line[2],
                                            current_line[3],
                                            str(0),
                                            current_line[5]])
                    virtual_bed_file = pybedtools.BedTool(bed_string, from_string=True)
                    result = exons.intersect(virtual_bed_file, s=True)
                    fasta_bed_line_start = ""
                    fasta_bed_line_stop = ""

                    start = 0
                    stop = 0

                    #print("Current line", current_line)
                    flanking_exon_cache[name] = {}
                    all_exons_circle[name] = []
                    for result_line in str(result).splitlines():
                        bed_feature = result_line.split('\t')
                        # this is a single-exon circRNA
                        #print("bed feature", bed_feature)
                        if bed_feature[1] == current_line[1] and bed_feature[2] == current_line[2]:
                            #print("Single exon condition")
                            # remove 1 bp from start and end to correct for the gap 
                            temp_bed_feature = bed_feature 
                            temp_bed_feature[1] = str(int(temp_bed_feature[1]) - 1)
                            #temp_bed_feature[2] = str(int(temp_bed_feature[2]) + 1)
                            result_line = "\t".join(temp_bed_feature)
                            #print("Updated result line", result_line)
                            fasta_bed_line_start += result_line + "\n"
                            start = 1
                            stop = 1
                            all_exons_circle[name].append([bed_feature[1], bed_feature[2]])

                        if bed_feature[1] == current_line[1] and start == 0:
                            #print("Start zero condition")
                            temp_bed_feature = bed_feature
                            temp_bed_feature[1] = str(int(temp_bed_feature[1]) - 1)
                            result_line = "\t".join(temp_bed_feature) 
                            fasta_bed_line_start += result_line + "\n"
                            start = 1
                            all_exons_circle[name].append([bed_feature[1], bed_feature[2]])

                        if bed_feature[2] == current_line[2] and stop == 0:
                            #print("Stop zero condition")
                            temp_bed_feature = bed_feature
                            temp_bed_feature[1] = str(int(temp_bed_feature[1]) - 1)
                            result_line = "\t".join(temp_bed_feature) 
                            fasta_bed_line_stop += result_line + "\n"
                            stop = 1
                            all_exons_circle[name].append([bed_feature[1], bed_feature[2]])

                        # these exons are kept for correctly drawing the circRNAs later
                        # not used for primer design
                        if bed_feature[1] > current_line[1] and bed_feature[2] < current_line[2]:
                            flanking_exon_cache[name][bed_feature[1] + "_" + bed_feature[2]] = 1
                            all_exons_circle[name].append([bed_feature[1], bed_feature[2]])

                    print(all_exons_circle)
                    print(name, all_exons_circle[name])
                    # first and last exons
                    virtual_bed_file_start = pybedtools.BedTool(fasta_bed_line_start, from_string=True)
                    virtual_bed_file_stop = pybedtools.BedTool(fasta_bed_line_stop, from_string=True)

                    virtual_bed_file_start = virtual_bed_file_start.sequence(fi=self.fasta_file, s=True)
                    virtual_bed_file_stop = virtual_bed_file_stop.sequence(fi=self.fasta_file, s=True)

                    if stop == 0 or start == 0:
                        print("Could not identify the exact exon-border of the circRNA.")
                        print("Will continue with non-annotated, manually extracted sequence.")
                        # we have to manually reset the start position

                        fasta_bed_line = "\t".join([current_line[0],
                                                    current_line[1],
                                                    current_line[2],
                                                    current_line[5]])
                        
                        virtual_bed_file_start = pybedtools.BedTool(fasta_bed_line, from_string=True)
                        virtual_bed_file_start = virtual_bed_file_start.sequence(fi=self.fasta_file, s=True)
                        virtual_bed_file_stop = ""
                    exon1 = ""
                    exon2 = ""

                    if virtual_bed_file_start:
                        exon1 = open(virtual_bed_file_start.seqfn).read().split("\n", 1)[1].rstrip()

                    if virtual_bed_file_stop:
                        exon2 = open(virtual_bed_file_stop.seqfn).read().split("\n", 1)[1].rstrip()

                    circ_rna_number += 1
                    print("extracting flanking exons for circRNA #", circ_rna_number, name, end="\n", flush=True)

                    # fetch the information about first/last circle that contributes to the BSJ
                    if current_line[5] == "+":
                        bsj_exon = all_exons_circle[name][-1]
                    elif current_line[5] == "-":
                        bsj_exon = all_exons_circle[name][0]
                    else:
                        print("No strand information present, assuming + strand")
                        bsj_exon = all_exons_circle[name][-1]
     
        else:
            print("Please provide Circtools detect output Coordinate file via option -d.")
            sys.exit(-1)
        
            if not exon_cache:
                print("Could not find any circRNAs matching your criteria, exiting.")
                exit(-1)
        print("Cleaning up")
        """      
        ## cleanup / delete tmp files
        os.remove(exon_storage_tmp)
        os.remove(blast_storage_tmp)
        os.remove(blast_xml_tmp)
        """