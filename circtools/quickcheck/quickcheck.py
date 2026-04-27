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
import subprocess
import circtools




class QuickCheck(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version
        self.command = sys.executable

    def module_name(self):
        """"Return a string representing the name of the module."""
        return self.program_name

    def run_module(self):

        # we first need to make sure all parameters are sane, the R script does not validate any parameters

        # check output directory

        import os

        # let's first check if the temporary directory exists
        if not os.path.exists(self.cli_params.output_directory):
            os.makedirs(self.cli_params.output_directory)

        # let's first check if the output directory exists
        if not (os.access(self.cli_params.output_directory, os.W_OK)):
            self.log_entry("Output directory %s not writable." % self.cli_params.output_directory)
            # exit with -1 error if we can't use it
            exit(-1)

        # check circtools detect directory
        if not (os.path.exists(self.cli_params.detect_dir)):
            self.log_entry("circtools detect data directory %s does not exist (or is not accessible)."
                           % self.cli_params.detect_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        # check circtools detect files (only existence, not the content)
        self.check_input_files([
                                self.cli_params.detect_dir+"CircRNACount",
                                self.cli_params.detect_dir+"LinearCount",
                                self.cli_params.detect_dir+"CircCoordinates"
                                ])

        # check STAR directory
        if not (os.path.exists(self.cli_params.star_dir)):
            self.log_entry("STAR data directory %s does not exist (or is not accessible)."
                           % self.cli_params.star_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        # check sample names
        if len(self.cli_params.condition_list.split(",")) < 2:
            self.log_entry("Error: Length of parameter list specified via -c is < 2.")

        # check grouping
        if len(self.cli_params.grouping.split(",")) < 2:
            self.log_entry("Error: Length of parameter list specified via -g is < 2.")

        for column in self.cli_params.grouping.split(","):
            try:
                int(column)
            except ValueError:
                self.log_entry("Error: group %s is no valid group index." % str(column))
                exit(-1)

        self.check_int_arguments([self.cli_params.remove_suffix_chars,
                                  self.cli_params.remove_prefix_chars])

        # needed for Rscript decoupling
        import subprocess

        # import re module
        import re

        python_location = subprocess.check_output(['which', self.command], universal_newlines=True,
                                                  stderr=subprocess.STDOUT).split('\n')[0]
        python_version = subprocess.check_output([self.command, '--version'], universal_newlines=True,
                                                 stderr=subprocess.STDOUT).strip()
        self.log_entry("Using Python interpreter %s [%s]" % (python_version, python_location))

        # ------------------------------------ need to call the correct R script here -----------------------

        # need to define path top R wrapper
        quickcheck_script = os.path.join(os.path.dirname(circtools.__file__), "scripts", "circtools_quickcheck_wrapper.py")



        if self.cli_params.starfolder == "":
            self.cli_params.starfolder = 0
            
        if not hasattr(self.cli_params, "remove_columns") or not self.cli_params.remove_columns:
            self.cli_params.remove_columns = "0"

        # Variable number of args in a list
        args = [
            self.cli_params.detect_dir,
            self.cli_params.star_dir,
            os.path.join(self.cli_params.output_directory, self.cli_params.output_name),
            self.cli_params.condition_list,
            self.cli_params.grouping,
            self.cli_params.colour,
            self.cli_params.cleanup,
            str(self.cli_params.starfolder),
            str(self.cli_params.remove_suffix_chars),
            str(self.cli_params.remove_prefix_chars),
            str(self.cli_params.remove_columns),
        ]

        cmd = [self.command, quickcheck_script] + args

        self.log_entry(f"DEBUG running: {' '.join(cmd)}")

        subprocess.run(cmd, check=True)
