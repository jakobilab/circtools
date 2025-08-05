#!/usr/bin/env python3

import os
import circ_module.circ_template
from scripts.circtools_circtest_wrapper import run_circ_test_wrapper  # your translated Python wrapper
import logging
os.makedirs("logs", exist_ok=True)


logging.basicConfig(
    filename="logs/circtest.log",        # log file path
    level=logging.INFO,                  # or DEBUG for more verbosity
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


class CircTest(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version

    def module_name(self):
        return self.program_name

    def run_module(self):
        # --- Input & Output checks ---
        if not os.path.exists(self.cli_params.output_directory):
            os.makedirs(self.cli_params.output_directory)

        if not os.access(self.cli_params.output_directory, os.W_OK):
            self.log_entry(f"Output directory {self.cli_params.output_directory} not writable.")
            exit(-1)

        if not os.path.exists(self.cli_params.detect_dir):
            self.log_entry(f"DCC/detect data directory {self.cli_params.detect_dir} does not exist.")
            exit(-1)

        self.check_input_files([
            os.path.join(self.cli_params.detect_dir, "CircRNACount"),
            os.path.join(self.cli_params.detect_dir, "LinearCount"),
            os.path.join(self.cli_params.detect_dir, "CircCoordinates")
        ])

        if len(self.cli_params.condition_list.split(",")) < 2:
            self.log_entry("Error: Length of condition list (-c) < 2.")
        if len(self.cli_params.condition_columns.split(",")) < 2:
            self.log_entry("Error: Length of condition columns (-l) < 2.")
        if len(self.cli_params.grouping.split(",")) < 2:
            self.log_entry("Error: Length of grouping (-g) < 2.")

        for val in self.cli_params.condition_columns.split(","):
            if not val.isdigit():
                self.log_entry(f"Error: column '{val}' is not a valid index.")
                exit(-1)

        for val in self.cli_params.grouping.split(","):
            if not val.isdigit():
                self.log_entry(f"Error: group '{val}' is not a valid group index.")
                exit(-1)

        self.check_int_arguments([
            self.cli_params.num_replicates,
            self.cli_params.filter_sample,
            self.cli_params.filter_count,
            self.cli_params.max_plots
        ])

        self.check_float_arguments([
            self.cli_params.max_fdr,
            self.cli_params.percentage,
            self.cli_params.range
        ])

        if not (0 < self.cli_params.max_fdr < 1):
            self.log_entry("Error: FDR (-f) must be in (0, 1).")
            exit(-1)

        if not (0 < self.cli_params.percentage <= 1):
            self.log_entry("Error: Percentage (-p) must be in (0, 1].")
            exit(-1)

        # --- Prepare arguments for Python wrapper ---
        class Args:
            pass

        args = Args()
        args.dcc_data = self.cli_params.detect_dir
        args.replicates = int(self.cli_params.num_replicates)
        args.condition_list = self.cli_params.condition_list
        args.condition_columns = self.cli_params.condition_columns
        args.output_name = os.path.join(self.cli_params.output_directory, self.cli_params.output_name)
        args.max_fdr = float(self.cli_params.max_fdr)
        args.max_plots = int(self.cli_params.max_plots)
        args.filter_sample = int(self.cli_params.filter_sample)
        args.filter_count = int(self.cli_params.filter_count)
        args.groups = self.cli_params.grouping
        args.output_label = self.cli_params.label
        args.percent_filter = float(self.cli_params.percentage)
        args.only_negative = bool(self.cli_params.only_negative)
        args.add_header = bool(self.cli_params.add_header)
        args.y_range = float(self.cli_params.range)
        args.colour_mode = self.cli_params.colour
        args.circle_description = min([int(c) for c in args.condition_columns.split(",")]) - 1


        # --- Run Python-native CircTest wrapper ---
        try:
            run_circ_test_wrapper(args)
        except Exception as e:
            self.log_entry(f"Error during CircTest execution: {e}")
            logging.error("Exception details:", exc_info=True)
            exit(-1)
