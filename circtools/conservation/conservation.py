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

# Loading the functionalities defined in other scripts as modules

# for now silence these warnings
from Bio import BiopythonDeprecationWarning
import warnings


with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from . import FetchOrthologs as FO
    from . import LiftOver as LO
    from . import FetchRegionGeneSequence as FS
    from . import SequenceAlignment as AL


import yaml
import circtools

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
        self.mm10_flag = self.cli_params.mm10_conversion
        self.hg19_flag = self.cli_params.hg19_conversion
        self.pairwise = self.cli_params.pairwise_flag
        # changes for config file addition
        self.config = self.cli_params.config

        # check if config present or is accessible
        if self.config:
            if os.path.isfile(self.config):
                print("Config file detected.")
            else:
                print("Configuration file not accessible. Using example config file.")
                self.config = os.path.dirname(os.path.realpath(__file__)) + "/config/example.config"
        else:
            print("Configuration argument not provided. Using example config file.")
            self.config = os.path.dirname(os.path.realpath(__file__)) + "/config/example.config"

        # process the available config file
        with open(self.config, 'r') as config_file:
            config = (yaml.safe_load(config_file))
            dict_species_ortholog = {}
            dict_species_liftover = {}
            dict_species_conservation = {}
            for each in config:
                # dictionary for fetchOrthologs.py script
                dict_species_ortholog[config[each]['ortho_id']] = config[each]['name']
                # dictionary for LiftOver.py script
                dict_species_liftover[config[each]['name']] = config[each]['id']
                # dictionary for conservation.py scipt
                dict_species_conservation[config[each]['input']] = config[each]['name']
            self.dict_species_ortholog = dict_species_ortholog
            self.dict_species_liftover = dict_species_liftover
            self.dict_species_conservation = dict_species_conservation

        # gene_list file argument
        if (self.cli_params.gene_list):
            self.gene_list = self.cli_params.gene_list
        elif (self.cli_params.gene_list_file):
            gene_file = [line.rstrip() for line in open(self.cli_params.gene_list_file[0])]
            self.gene_list = gene_file
        else:
            print("Need to provide gene list by either -G or -GL options!")
            exit(-1)

        # target species argument
        if (self.cli_params.target_species):
            # argument is comma separated list of target species for which conservation will be calculated
            self.target_species = self.cli_params.target_species.split(",")
            if any(e not in ["mm", "hs", "ss", "rn", "cl"] for e in self.target_species):
                print("Please only specify available choices of organisms: mm, hs, ss, rn or cl")
                exit(-1)
        else:
            print("Need to provide target species to calculate conservation")
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
        except (PermissionError, FileNotFoundError):
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

    def process_data(self, ):
        # reading information from config file
        with open(configuration, 'r') as config_file:
            config = (yaml.safe_load(config_file))
    
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
        letters = string.ascii_letters
        tmp_prefix =  ''.join(random.choice(letters) for i in range(10))        # prefix for tmp files
              
        # species_dictionary = {"mm": "mouse", "hs": "human", "rn": "rat", "ss": "pig", "cl": "dog"}
        species_dictionary = self.dict_species_conservation
        dict_species_ortholog = self.dict_species_ortholog
        dict_species_liftover = self.dict_species_liftover

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

                    # if mouse and human assemblies are told, convert the the co-ordinates using liftover function
                    if self.mm10_flag:
                        #print("Before liftover:", current_line)
                        print("Converting mm10 coordinates to mm39 assembly")
                        lifted = LO.liftover("mouse", "mouse", current_line, self.temp_dir, tmp_prefix, {}, "mm10", self.dict_species_liftover)
                        current_line = lifted.parseLiftover()
                    elif self.hg19_flag:
                        #print("Before liftover:", current_line)
                        print("Converting hg19 coordinates to hg38 assembly") 
                        lifted = LO.liftover("human", "human", current_line, self.temp_dir, tmp_prefix, {}, "hg19", self.dict_species_liftover)
                        current_line = lifted.parseLiftover()

                    sep = "_"
                    name = sep.join([current_line[3],
                                        current_line[0],
                                        current_line[1],
                                        current_line[2],
                                        current_line[5]])
                    print(name)
                    if self.id_list and not self.gene_list and name not in self.id_list:
                        continue

                    #circrna_length = int(current_line[2]) - int(current_line[1])

                    # check for intronic circles.
                    if (current_line[6] == "intron-intron"):
                        print("Warning!!! This is an intronic circle. The conservation will not be checked.")
                        continue

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
                        if bed_feature[1] == current_line[1] and bed_feature[2] == current_line[2]:
                            # remove 1 bp from start and end to correct for the gap 
                            temp_bed_feature = bed_feature 
                            temp_bed_feature[1] = str(int(temp_bed_feature[1]) - 1)
                            result_line = "\t".join(temp_bed_feature)
                            fasta_bed_line_start += result_line + "\n"
                            start = 1
                            stop = 1
                            all_exons_circle[name].append([bed_feature[1], bed_feature[2]])

                        if bed_feature[1] == current_line[1] and start == 0:
                            temp_bed_feature = bed_feature
                            temp_bed_feature[1] = str(int(temp_bed_feature[1]) - 1)
                            result_line = "\t".join(temp_bed_feature) 
                            fasta_bed_line_start += result_line + "\n"
                            start = 1
                            all_exons_circle[name].append([bed_feature[1], bed_feature[2]])

                        if bed_feature[2] == current_line[2] and stop == 0:
                            temp_bed_feature = bed_feature
                            temp_bed_feature[1] = str(int(temp_bed_feature[1]) - 1)
                            result_line = "\t".join(temp_bed_feature) 
                            fasta_bed_line_stop += result_line + "\n"
                            stop = 1
                            all_exons_circle[name].append([bed_feature[1], bed_feature[2]])

                    #print(all_exons_circle)

                    # first and last exons
                    virtual_bed_file_start = pybedtools.BedTool(fasta_bed_line_start, from_string=True)
                    virtual_bed_file_stop = pybedtools.BedTool(fasta_bed_line_stop, from_string=True)

                    virtual_bed_file_start = virtual_bed_file_start.sequence(fi=self.fasta_file, s=True)
                    virtual_bed_file_stop = virtual_bed_file_stop.sequence(fi=self.fasta_file, s=True)

                    # print("Virtual bed start and stop", virtual_bed_file_start, virtual_bed_file_stop)

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


                    fasta_out_string = ""               # the fasta file string for sequennce alignment
                    circ_sequence = exon2 + exon1       # this is the joint exon circular RNA sequence to use for alignment
                    fasta_out_string = fasta_out_string + ">" + self.organism + "(" + current_line[0] + ":" + current_line[1] + "-" + current_line[2] + ")" + "\n" + circ_sequence + "\n"  

                    #print("All exons circle: ", all_exons_circle[name])
                    # fetch the information about first/last circle that contributes to the BSJ
                    print("extracting flanking exons for circRNA #", circ_rna_number, name, end="\n", flush=True)

                    if len(all_exons_circle[name]) == 1:
                        # this is a single exon circle
                        bsj_exon = all_exons_circle[name][0] 
                    else:
                        if current_line[5] == "+":
                            bsj_exon = all_exons_circle[name][-1]
                            first_exon = all_exons_circle[name][0]
                        elif current_line[5] == "-":
                            print(all_exons_circle[name])
                            bsj_exon = all_exons_circle[name][0]
                            first_exon = all_exons_circle[name][-1]
                        else:
                            print("No strand information present, assuming + strand")
                            bsj_exon = all_exons_circle[name][-1]
                            first_exon = all_exons_circle[name][0]

                    # call exon fetchorthologs function to store orthologs which then will be sent to liftover function
                    host_gene = current_line[3]
                    fetchOrtho = FO.fetch(host_gene, species_dictionary[self.organism], self.dict_species_ortholog)
                    ortho_dict = fetchOrtho.parse_json()
                    # check if each target species is present in the dictionary and remove if not
                    for each_species in self.target_species:
                        species_name = species_dictionary[each_species]
                        if species_name not in ortho_dict.keys():
                            print("No ortholog found in species" + species_name + ". Removing this species for further analysis.")
                            self.target_species.remove(each_species)
                    # after removing this species, check if there are any species left to perform the analysis
                    if len(self.target_species) == 0:
                        print("No species found to perform conservation analysis.")
                        sys.exit()

                    #print(current_line)

                    for each_target_species in self.target_species:
                        print("Processing target species: ",each_target_species)             
                    
                        # take these flanking exons per circle and perform liftover 
                        if "first_exon" in locals():
                            # liftover first exon
                            print("*** Lifting over first exon ***")
                            first_line = [current_line[0], first_exon[0], first_exon[1], current_line[3], current_line[4], current_line[5]]
                            lifted = LO.liftover(species_dictionary[self.organism], species_dictionary[each_target_species], first_line, self.temp_dir, tmp_prefix+"_first", ortho_dict, "other", self.dict_species_liftover)
                            first_exon_liftover = lifted.find_lifted_exons()
                            if (first_exon_liftover == None):
                                print("No lifted co-ordinates found for first exon. Skipping " + species_dictionary[each_target_species] + "for further analysis.")
                                continue
                        
                        print("*** Lifting over BSJ exon ***")
                        bsj_line = [current_line[0], bsj_exon[0], bsj_exon[1], current_line[3], current_line[4], current_line[5]]
                        lifted = LO.liftover(species_dictionary[self.organism], species_dictionary[each_target_species], bsj_line, self.temp_dir, tmp_prefix+"_BSJ", ortho_dict, "other", self.dict_species_liftover)
                        bsj_exon_liftover = lifted.find_lifted_exons()
                        if (bsj_exon_liftover == None):
                            print("No lifted co-ordinates found for BSJ exon. Skipping " + species_dictionary[each_target_species] + "for further analysis.")
                            continue


                        # fetch sequences for both these exons now
                        if "first_exon" in locals():
                            first_exon_seq = FS.sequence(species_dictionary[each_target_species], first_exon_liftover)
                            bsj_exon_seq = FS.sequence(species_dictionary[each_target_species], bsj_exon_liftover)
                            #print("Lifted over coordinates:", first_exon_liftover, bsj_exon_liftover)
                            lifted_circle = first_exon_liftover[:2] + [bsj_exon_liftover[2]]
                            if current_line[5] == "+":
                                circ_sequence_target = str(bsj_exon_seq.fetch_sequence()) + str(first_exon_seq.fetch_sequence())
                            elif current_line[5] == "-":
                                circ_sequence_target = str(first_exon_seq.fetch_sequence()) + str(bsj_exon_seq.fetch_sequence())
                            else:
                                circ_sequence_target = str(bsj_exon_seq.fetch_sequence()) + str(first_exon_seq.fetch_sequence())
                        else:
                            bsj_exon_seq = FS.sequence(species_dictionary[each_target_species], bsj_exon_liftover)
                            circ_sequence_target = str(bsj_exon_seq.fetch_sequence())
                            #print("Lifted over coordinates:", bsj_exon_liftover)
                            lifted_circle = bsj_exon_liftover

                        lifted_circle = list(map(str, lifted_circle))
                        print("Lifted circle in target species ", each_target_species , " is " , lifted_circle)
                        out_bed_file = self.output_dir + "lifted_" + name + "_" +   each_target_species + ".bed"
                        with open(out_bed_file, "w") as bed_out:
                            bed_out.write("\t".join(lifted_circle)+"\n")        

                        # writing into fasta file for alignments
                        fasta_out_string = fasta_out_string + ">" + each_target_species + "(" + lifted_circle[0] + ":" + lifted_circle[1] + "-" + lifted_circle[2] + ")" + "\n" + circ_sequence_target + "\n" 
                        
                    fasta_file_alignment = self.output_dir + "/alignment_" + name + ".fasta"
                    with open(fasta_file_alignment, "w") as fasta_out:
                        fasta_out.write(fasta_out_string)

                    # call multiple sequencce alignment function
                    align = AL.Alignment(fasta_file_alignment, self.organism, name, self.output_dir)
                    align.draw_phylo_tree()
                    # if pairwise alignment flag is turned on, run the pairwise alignment function
                    if (self.pairwise):
                        align.pairwise_alignment()

        else:
            print("Please provide Circtools detect output Coordinate file via option -d.")
            sys.exit(-1)
        
        print("Cleaning up")
        """      
        ## cleanup / delete tmp files
        os.remove(exon_storage_tmp)
        os.remove(blast_storage_tmp)
        os.remove(blast_xml_tmp)
        """
        