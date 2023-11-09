#!/usr/bin/env python3
# @Author: Shubhada Kulkarni
# @Email:  shubhada.kulkarni@uni-heidelberg.de
# @Institute: University Hospital Heidelberg, Heidelberg


# time python combine_circtools_ciriquant.py -c /scratch/circtools2/circtools/sample_data/CircRNACount -l /scratch/circtools2/circtools/sample_data/LinearCount -r /scratch/circtools2/circtools/sample_data/CircCoordinates -q list_ciriquant_sampledata -o test -s __R1.Chimeric.out.junction

##########################################################################
# IMPORTING REQUIRED LIBRARIES
import re
import csv
import sys
import os
import glob
import warnings
import itertools
from math import ceil
from optparse import OptionParser
from operator import itemgetter
import pandas as pd
import numpy as np
import itertools
import logging

##########################################################################
class metatool():

    def __init__(self):
        #self.tmp_dir = tmp_dir
        print("Merging counts from Circtools and CIRIquant\n")

    '''
    # Commenting because this part will be checked in the detect module class
    @staticmethod
    def parseoptions():
        # COMMAND LINE ARGUMENTS
        parser=OptionParser()
        parser.add_option('-c', "--circtools", help="CircRNACount output file of Circtools detect step")
        parser.add_option('-l', "--linear", help="LinearCount output file of Circtools detect step")
        parser.add_option('-r', "--circcoordinates", help="CircCoordinate output file of Circtools detect step. Required to fetch strand information")
        parser.add_option('-q', "--ciriquant", help="File with CIRIquant output files per sample. One sample output file path on one line, ordered as columns in CircRNACount file")
        parser.add_option('-s', "--string", help="ID/String to add in the output file name")
        parser.add_option('-o', "--out_dir", help="Output dictionary name") #, default="combined_outputs_Circtools-CIRIquant")
        (opts, args)=parser.parse_args()
        return(opts)
    
    @staticmethod
    def check_python_version():
       # CHECKING PYTHON VERSION
        MIN_PYTHON = (3, 0)
        if sys.version_info < MIN_PYTHON:
            sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

    def check_command_options(self, opts):
        # OPTIONS CHECK
        if (opts.circtools == None):
            print("Please provide Circtools Detect output file CircRNACount file using option -c")
            exit(1)

        if (opts.linear == None):
            print("Please provide Circtools Detect output file LinearCount file using option -l")
            exit(1)

        if (opts.circcoordinates == None):
            print("Please provide Circtools Detect output file CircCoordinates file using option -co")
            exit(1)

        if (opts.ciriquant == None):
            print("Please provide list of CIRIquant output in one file using option -r")
            exit(1)

        if (opts.out_dir == None):
            print("Please provide output directory name using option -o")
            exit(1)
        else:
            if (os.path.isdir(opts.out_dir)):
                print("Reusing existing directory: " + opts.out_dir)
            else:
                os.mkdir(opts.out_dir)
        return(opts)
    '''
    ##########################################################################
    # INPUT FILE READING
    # Read and store circtools output

    def parseCirctools(self, input_circtools, input_circcordinates, input_linear, replace_string):
        # function to parse the output of circtools and store it in a dictionary
        self.input_circtools = input_circtools
        self.input_circcordinates = input_circcordinates
        self.input_linear = input_linear
        self.replace_string = replace_string

        dict_circtools = {}
        dict_circtools_linear = {}
        inp = open(input_circtools, mode='r').readlines()
        inp_co = open(input_circcordinates, mode='r').readlines()
        inp_l = open(input_linear, mode='r').readlines()
        num_c = len(inp)
        num_co = len(inp_co)
        num_l = len(inp_co)
        if (num_c != num_co):
            print("Number of lines in CircRNACount and CircCoordinates are different!")
            exit(1)
        if (num_c != num_l):
            print("Number of lines in CircRNACount and LinearCount are different!")
            exit(1)
        header = inp[0].rstrip().split('\t')
        samplenames_circtools =  [s.replace(replace_string, "") for s in header[3:]]
        #print(samplenames_circtools)
        index = 1
        for line in inp[1:]:
            line = line.rstrip().split('\t')
            circ_id = line[0] + "_" + line[1] + "_" + line[2]

            line_l = inp_l[index].strip().split('\t')
            circ_id_l = line_l[0] + "_" + line_l[1] + "_" + line_l[2]

            # fetch strand info from CircCoordinates file
            line_2 = inp_co[index].strip().split('\t')
            circ_id_2 = line_2[0] + "_" + line_2[1] + "_" + line_2[2]
            if (circ_id != circ_id_2):
                print("Lines in CircRNACount and CircCoordinates do not match!")
                print(line)
                print(line_2)
                exit(1)
            elif (circ_id != circ_id_l):
                print("Lines in CircRNACount and LinearCount do not match!")
                print(line)
                print(line_l)
                exit(1)
            else:
                # update circ_id with strand information
                circ_id = circ_id + "_" + line_2[5]
            dict_circtools[circ_id] = list(map(int, line[3:]))
            dict_circtools_linear[circ_id] = list(map(int, line_l[3:]))
            
            index = index + 1

        ## Normalization of circtools circular and linear dataframes
        pd_circtools = pd.DataFrame.transpose(pd.DataFrame.from_dict(dict_circtools))                         # non-normalized circtools circular counts
        pd_circtools_linear = pd.DataFrame.transpose(pd.DataFrame.from_dict(dict_circtools_linear))                  # non-normalized circtools linear count
        pd_circtools.columns = samplenames_circtools                                                                 
        pd_circtools_linear.columns = samplenames_circtools                                                                 
        normalized_pd_circtools=(pd_circtools-pd_circtools.min())/(pd_circtools.max()-pd_circtools.min())       # normalized circtools circular counts
        normalized_pd_circtools_linear=(pd_circtools_linear-pd_circtools_linear.min())/(pd_circtools_linear.max()-pd_circtools_linear.min())       # normalized circtools linear counts
        #print(normalized_pd_circtools_linear)
        return (pd_circtools, pd_circtools_linear, normalized_pd_circtools, normalized_pd_circtools_linear)

    ##########################################################################
    # Create a two-dimentional dictionary of CIRIquant matches. 
    def parseCIRIquant(self, list_ciriquant):
        # function to parse CIRIquant output files one-by-one for all samples in the input file
        self.list_ciriquant = list_ciriquant
        dict_ciriquant = {}
        dict_ciriquant_linear = {}
        dict_ciriquant_norm = {}
        dict_ciriquant_linear_norm = {}
        inp = open(list_ciriquant, "r").readlines()
        for line in inp:
                line = line.rstrip().split("\t")
                #print(line)
                if (os.path.isfile(line[1]) == False):
                    print("Following CIRIquant output file does not exist: " + line[1])
                    exit(1)
                
                samplename = line[0]
                #print(samplename)
                fin = open(line[1]).readlines()
                dict_ciriquant[samplename] = {}
                dict_ciriquant_linear[samplename] = {}
                for eachline in fin:
                    if ("circRNA_ID" in eachline):
                        continue
                    else:
                        eachline = eachline.rstrip().split("\t")
                    key_string = eachline[1] + "_" + eachline[2] + "_" + eachline[3] + "_" + eachline[10]
                    dict_ciriquant[samplename][key_string] = [int(eachline[4])]
                    doptions.list_ciriquant, output_circ_counts, output_coordinates, output_linear, 
                    options.cleanup, options.out_dirict_ciriquant_linear[samplename][key_string] = [int(eachline[6])]

                cq_pd = pd.DataFrame.transpose(pd.DataFrame.from_dict(dict_ciriquant[samplename]))
                cq_pd_linear = pd.DataFrame.transpose(pd.DataFrame.from_dict(dict_ciriquant_linear[samplename]))
                dict_ciriquant[samplename] = cq_pd
                dict_ciriquant_linear[samplename] = cq_pd_linear
                normalized_pd_ciriquant=(cq_pd)/(cq_pd.max())
                normalized_pd_ciriquant_linear=(cq_pd_linear)/(cq_pd_linear.max())
                dict_ciriquant_norm[samplename] = normalized_pd_ciriquant
                dict_ciriquant_linear_norm[samplename] = normalized_pd_ciriquant_linear
        return (dict_ciriquant, dict_ciriquant_linear, dict_ciriquant_norm, dict_ciriquant_linear_norm)

    ############################################################################
    # MERGING AND WRITING OUTPUTS TO A FILE
    def merging(self, ciriquant, circtools, circcoordinates, linear, string, out_dir):

        # parse the outputs of both tools
        out_ciriquant = self.parseCIRIquant(ciriquant)
        out_circtools = self.parseCirctools(circtools, circcoordinates, linear, string)

        # Output dataframes
        temp_cq = [x.index.values for x in out_ciriquant[0].values()]
        all_matches = list(set(list(set(out_circtools[0].index.values)) + list(set([j for i in temp_cq for j in i]))))

        combined_circ = pd.DataFrame(0, index=all_matches, columns=out_circtools[0].columns)
        combined_linear = pd.DataFrame(0, index=all_matches, columns=out_circtools[0].columns)
        combined_circ_norm = pd.DataFrame(0.0, index=all_matches, columns=out_circtools[0].columns)
        combined_linear_norm = pd.DataFrame(0.0, index=all_matches, columns=out_circtools[0].columns)

        # Go over every element and fill in the values from circtools and CIRIquant. 
        for i in combined_circ.index.values:
            for j in combined_circ.columns.tolist():
                present_ctools = 0
                present_cquant = 0
                # first, check if the circle is predicted by both tools
                try:
                    score_ctools = out_circtools[0].at[i,j]
                    present_ctools = 1
                except:
                    pass

                try:
                    score_cquant = out_ciriquant[0][j].at[i,0]
                    present_cquant = 1
                except:
                    pass

                # if the circle is predicted by both, print the details for later analysis
                if ((present_ctools == 1) and (present_cquant == 1)):
                    combined_circ.at[i,j] = out_circtools[0].at[i,j]
                    combined_linear.at[i,j] = out_circtools[1].at[i,j]
                    combined_circ_norm.at[i,j] = out_circtools[2].at[i,j]
                    combined_linear_norm.at[i,j] = out_circtools[3].at[i,j]
                elif (present_ctools == 1):
                    combined_circ.at[i,j] = out_circtools[0].at[i,j]
                    combined_linear.at[i,j] = out_circtools[1].at[i,j]
                    combined_circ_norm.at[i,j] = out_circtools[2].at[i,j]
                    combined_linear_norm.at[i,j] = out_circtools[3].at[i,j]
                elif (present_cquant == 1):
                    combined_circ.at[i,j] = out_ciriquant[0][j].at[i,0]
                    combined_linear.at[i,j] = out_ciriquant[1][j].at[i,0]
                    combined_circ_norm.at[i,j] = out_ciriquant[2][j].at[i,0]
                    combined_linear_norm.at[i,j] = out_ciriquant[3][j].at[i,0]


        # WRITING OUTPUTS TO A FILES
        combined_circ.to_csv(out_dir + "/CircRNACount_Merged", sep = "\t")
        combined_linear.to_csv(out_dir + "/LinearCount_Merged", sep = "\t")
        combined_circ_norm.to_csv(out_dir + "/CircRNACount_Merged_Normalized", sep = "\t")
        combined_linear_norm.to_csv(out_dir + "/LinearCount_Merged_Normalized", sep = "\t")

        print("Circtools and CIRIquant matches merged successfully")
