#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -J "circtools reconstruct"
#SBATCH -p long
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Sample name] [target dir e.g. /path/to/data/] [BED file] [DCC dirrectory] [CircRNACount directory]"
  exit
fi

#module load circtools
#module load samtools

main_out=$2/
sample_name=$1
bed_file=$3
dcc_dir=$4
tmp_folder=/scratch/global_tmp/
dcc_out_dir=$5

#######################

main_bam=$dcc_dir/${sample_name}.bam
main_junction=$dcc_dir/${sample_name}.Chimeric.out.junction

mate1_bam=$dcc_dir/${sample_name}.mate1.bam
mate1_junction=$dcc_dir/${sample_name}.mate1.Chimeric.out.junction

mate2_bam=$dcc_dir/${sample_name}.mate2.bam
mate2_junction=$dcc_dir/${sample_name}.mate2.Chimeric.out.junction.fixed

merged_bam=$main_out/${sample_name}_merged.bam

###################

# preprocessing

# merge both mate BAM files into one new BAM file
samtools merge -l 9 -@ 40 $merged_bam $main_bam $mate1_bam $mate2_bam

# re-index the newly aggregated BAM file
samtools index $merged_bam

FUCHS -N $sample_name -D $dcc_out_dir/CircRNACount -B $merged_bam -A $bed_file -O $main_out -F $mate2_junction -R $mate2_junction -J $main_junction -T $tmp_folder -p ensembl -r 2 -e 3 -q 2 -P 40

