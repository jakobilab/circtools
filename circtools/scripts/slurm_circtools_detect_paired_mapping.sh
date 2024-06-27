#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -J "circtools alignment"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

#module load star

# check if we have 6 arguments
if [ ! $# == 6 ]; then
  echo "Usage: $0 [STAR index] [Read 1 file] [Read 2 file] [target dir e.g. /tmp/] [Read 1 marker, e.g. R1] [GTF file]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2/$5/} : '\(.*\)\..*\.'`

# create the target directory, STAR will not do that for us
mkdir -pv $4/$target
mkdir -pv $4/$target/mate1/
mkdir -pv $4/$target/mate2/

# create random string
TMP_RND=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 10 | head -n 1`

#mkdir -pv /scratch/global_tmp/${TMP_RND}_${target}
#mkdir -pv /scratch/global_tmp/${TMP_RND}_${target}_mate1
#mkdir -pv /scratch/global_tmp/${TMP_RND}_${target}_mate2

OLD_PATH=`pwd`

# main mapping part

STAR	--runThreadN 40\
	--genomeDir $1\
	--genomeLoad NoSharedMemory\
	--outTmpDir /scratch/global_tmp/${TMP_RND}_${target}/\
	--readFilesIn $2 $3\
	--readFilesCommand zcat\
	--outFileNamePrefix $4/$target/\
	--outReadsUnmapped Fastx\
	--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
	--outSJfilterOverhangMin 15   15   15   15\
	--outFilterMultimapNmax 20\
	--chimMultimapNmax 20\
	--outFilterScoreMin 1\
	--outFilterMatchNminOverLread 0.7\
	--outFilterMismatchNmax 999\
	--outFilterMismatchNoverLmax 0.05\
	--alignIntronMin 20\
	--alignIntronMax 1000000\
	--alignMatesGapMax 1000000\
	--alignSJoverhangMin 15\
	--alignSJDBoverhangMin 10\
	--alignSoftClipAtReferenceEnds No\
	--chimSegmentMin 15\
	--chimScoreMin 15\
	--chimScoreSeparation 10\
	--chimJunctionOverhangMin 15\
	--sjdbGTFfile $6\
	--quantMode GeneCounts\
	--twopassMode Basic\
	--chimOutType Junctions

cd $4/$target

gzip Unmapped.out.mate1
gzip Unmapped.out.mate2

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam

grep "^@" Aligned.out.sam > header.txt

rm -f Aligned.out.sam
rm -f Chimeric.out.sam

rm -f -r _STARgenome
rm -f -r _STARpass1

samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
mv Aligned.noS.tmp Aligned.noS.bam
samtools index Aligned.noS.bam

samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin
samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp
mv Chimeric.noS.tmp Chimeric.noS.bam
samtools index Chimeric.noS.bam

rm -f Aligned.noS.sam
rm -f Chimeric.noS.sam

cd $OLD_PATH

## done with main mapping

## mapping mate1 now


STAR	--runThreadN 40\
	--genomeDir $1\
	--genomeLoad NoSharedMemory\
	--outTmpDir /scratch/global_tmp/${TMP_RND}_${target}_mate1/\
	--readFilesIn $4/$target/Unmapped.out.mate1.gz\
	--readFilesCommand zcat\
	--outFileNamePrefix $4/$target/mate1/ \
	--outReadsUnmapped Fastx\
	--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
	--outSJfilterOverhangMin 15   15   15   15\
	--outFilterMultimapNmax 20\
	--chimMultimapNmax 20\
	--outFilterScoreMin 1\
	--outFilterMatchNminOverLread 0.7\
	--outFilterMismatchNmax 999\
	--outFilterMismatchNoverLmax 0.05\
	--alignIntronMin 20\
	--alignIntronMax 1000000\
	--alignMatesGapMax 1000000\
	--alignSJoverhangMin 15\
	--alignSJDBoverhangMin 10\
	--alignSoftClipAtReferenceEnds No\
	--chimSegmentMin 15\
	--chimScoreMin 15\
	--chimScoreSeparation 10\
	--chimJunctionOverhangMin 15\
	--sjdbGTFfile $6\
	--quantMode GeneCounts\
	--twopassMode Basic\
	--chimOutType Junctions

cd $4/$target/mate1/

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam

grep "^@" Aligned.out.sam > header.txt

rm -f Aligned.out.sam
rm -f Chimeric.out.sam

rm -f -r _STARgenome
rm -f -r _STARpass1

samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
mv Aligned.noS.tmp Aligned.noS.bam
samtools index Aligned.noS.bam

samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin
samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp
mv Chimeric.noS.tmp Chimeric.noS.bam
samtools index Chimeric.noS.bam

rm -f Aligned.noS.sam
rm -f Chimeric.noS.sam

cd $OLD_PATH

## done with mate1 mapping

## mapping mate2 now

STAR	--runThreadN 40\
	--genomeDir $1\
	--genomeLoad NoSharedMemory\
	--outTmpDir /scratch/global_tmp/${TMP_RND}_${target}_mate2/\
	--readFilesIn $4/$target/Unmapped.out.mate2.gz\
	--readFilesCommand zcat\
	--outFileNamePrefix $4/$target/mate2/ \
	--outReadsUnmapped Fastx\
	--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
	--outSJfilterOverhangMin 15   15   15   15\
	--outFilterMultimapNmax 20\
	--chimMultimapNmax 20\
	--outFilterScoreMin 1\
	--outFilterMatchNminOverLread 0.7\
	--outFilterMismatchNmax 999\
	--outFilterMismatchNoverLmax 0.05\
	--alignIntronMin 20\
	--alignIntronMax 1000000\
	--alignMatesGapMax 1000000\
	--alignSJoverhangMin 15\
	--alignSJDBoverhangMin 10\
	--alignSoftClipAtReferenceEnds No\
	--chimSegmentMin 15\
	--chimScoreMin 15\
	--chimScoreSeparation 10\
	--chimJunctionOverhangMin 15\
	--sjdbGTFfile $6\
	--quantMode GeneCounts\
	--twopassMode Basic\
	--chimOutType Junctions

cd $4/$target/mate2/

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam

grep "^@" Aligned.out.sam > header.txt

rm -f Aligned.out.sam
rm -f Chimeric.out.sam

rm -f -r _STARgenome
rm -f -r _STARpass1

samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
mv Aligned.noS.tmp Aligned.noS.bam
samtools index Aligned.noS.bam

samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin
samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp
mv Chimeric.noS.tmp Chimeric.noS.bam
samtools index Chimeric.noS.bam

rm -f Aligned.noS.sam
rm -f Chimeric.noS.sam

# remove tmp dirs

rm /scratch/global_tmp/${TMP_RND}_${target}/ -rf
rm /scratch/global_tmp/${TMP_RND}_${target}_mate1/ -rf
rm /scratch/global_tmp/${TMP_RND}_${target}_mate2/ -rf

