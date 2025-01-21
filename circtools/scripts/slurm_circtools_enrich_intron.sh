#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 30
#SBATCH --mem=200G
#SBATCH -J "circtools enrich"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

#module load circtools

# check if we have 8 arguments
if [ ! $# == 8 ]; then
  echo "Usage: $0 [Chromosome sizes file] [Genome GTF file] [CLIP peak BED file] [CircRNA BED files] [output prefix] [output path] [# iterations] [TMP dir]"
  exit
fi

/usr/bin/time -v circtools enrich -c ${4} -b ${3} -a ${2} -g ${1} -i ${7} -I intron -p 20 -P 1 -T 1 -o ${6}/ -F ${5} -t ${8}/
