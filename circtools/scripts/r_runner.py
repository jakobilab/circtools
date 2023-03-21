#!/usr/bin/env python3

# Copyright (C) 2023 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# runner idea from https://stackoverflow.com/a/73649412/1900920

import subprocess
import sys
import circtools as _

base_path = _.__path__[0]


def _run(r_script):
    return subprocess.call(r_script, shell=True)


def circtools_circtest_wrapper():
    return _run("Rscript " +
                base_path + "/scripts/circtools_circtest_wrapper.R " +
                " ".join(sys.argv[1:]))


def circtools_enrich_visualization():
    return _run("Rscript " +
                base_path + "/scripts/circtools_enrich_visualization.R " +
                " ".join(sys.argv[1:]))


def circtools_exon_wrapper():
    return _run("Rscript " +
                base_path + "/scripts/circtools_exon_wrapper.R " +
                " ".join(sys.argv[1:]))


def circtools_primex_formatter():
    return _run("Rscript " +
                base_path + "/scripts/circtools_primex_formatter.R " +
                " ".join(sys.argv[1:]))


def circtools_primex_wrapper():
    return _run("Rscript " +
                base_path + "/scripts/circtools_primex_wrapper.R " +
                " ".join(sys.argv[1:]))


def circtools_quickcheck_wrapper():
    return _run("Rscript " +
                base_path + "/scripts/circtools_quickcheck_wrapper.R " +
                " ".join(sys.argv[1:]))


def circtools_reconstruct_visualization():
    return _run("Rscript " +
                base_path + "/scripts/circtools_reconstruct_visualization.R " +
                " ".join(sys.argv[1:]))


def circtools_sirna_formatter():
    return _run("Rscript " +
                base_path + "/scripts/circtools_sirna_formatter.R " +
                " ".join(sys.argv[1:]))

def circtools_reconstruct_summarized_coverage_profiles():
    return _run("Rscript " +
                base_path +
                "/scripts/circtools_reconstruct_summarized_coverage_profiles.R " +
                " ".join(sys.argv[1:]))

def circtools_reconstruct_coverage_graph():
    return _run("Rscript " +
                base_path + "/scripts/circtools_reconstruct_coverage_graph.R " +
                " ".join(sys.argv[1:]))


def install_R_dependencies():
    return _run("Rscript " + base_path
                + "/scripts/install_R_dependencies.R " + base_path)
