#!/bin/bash

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

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -J "Stringtie"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

# module load stringtie

# check if we have 3 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [BAM file] [GTF file] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> BAM file
# $2 -> GTF file
# $3 -> output directory

# create the target directory
mkdir $3 -p
#rm $3/*.gtf

stringtie $1 -v -f 0.2 -p 10 -G $2 -e -B -o $3/ballgown.gtf -C $3/reference_transcripts_full_coverage.gtf

#stringtie $1 -p 20 -G $2 -B -e 
#stringtie $1 -p 40 -G $2 -e -o $3/stringtie.gtf -C $3/reference_transcripts_full_coverage.gtf -B
#stringtie $1 -p 20 -G $2 -o $3/stringtie.gtf -j 2 -C $3/reference_transcripts_full_coverage.gtf
