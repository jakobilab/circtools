#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -J "bowtie2 rRNA filtering"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [rRNA index argument] [Read 1 file] [Read 2 file] [target dir e.g. /awesome/project/] [R1 marker, e.g. R1 or 1_sequence]"
  exit
fi

# $1 -> rRNA index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory
# $5 -> R1 marker

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2/$5/} : '\(.*\)\..*\.' | sed 's/_1$//g' `
#echo $target
#exit
# load the bowtie2 module
module load bowtie2

# SAM output goes to /dev/null
# run on 20 CPUs
# set fixed seed
# memory mapped IO for multiple instances
# display timing information
# write gz unmapping reads [== no rRNA] to target dir

bowtie2 -S /dev/null -x $1 -1 $2 -2 $3 --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-conc-gz $4/$target.fastq.gz 2> $4/$target.log
#bowtie2 -x $1 -1 $2 -2 $3 -S $4/$target.sam --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-conc-gz $4/$target.fastq.gz 2> $4/$target.log
