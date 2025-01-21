#!/bin/bash

# @Author: Shubhada Kulkarni <tjakobi>
# @Email:  shubhada.kulkarni@uni-heidelberg.de
# @Institute: University Hospital Heidelberg, Heidelberg

#SBATCH -J CIRIquant_run
#SBATCH --mail-type=END
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"
#SBATCH -n 8                # Number of cores
   
# -------------------------------

# Check if we have two arguments -> the input fastq files for paired-end
# sequencing
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Read 1 file] [Read 2 file] [target dir e.g. /awesome/project/] [config.yml file]"
  exit
fi

fastq1=$1     # Read 1 fastq file
fastq2=$2     # Read 2 fastq file
out=$3        # output directory
config=$4     # config file required for CIRIquant

echo $fastq1, $fastq2, $out, $config

CIRIquant -t 16 -1 $fastq1 -2 $fastq2 --config $config -o $out -p $out 
