#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -J "prepare circtools"

# check if we have 2 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [STAR source dir] [circtools destination dir]"
  exit
fi


SRC=${1}
DEST=${2}

if [ ! -d "$SRC" ]; then
  echo "Source directory $SRC does not exist!"
  exit;
fi

if [ ! -d "$DEST" ]; then
  echo "DCC directory $DEST does not exist!"
  exit;
fi

cd $SRC/

parallel ln -s `pwd`/{1}/mate{2}/Chimeric.out.junction ../$DEST/{1}.mate{2}.Chimeric.out.junction ::: * ::: 1 2
parallel ln -s `pwd`/{1}/mate{2}/Aligned.noS.bam ../$DEST/{1}.mate{2}.bam ::: * ::: 1 2 
parallel ln -s `pwd`/{1}/mate{2}/Aligned.noS.bam.bai ../$DEST/{1}.mate{2}.bam.bai ::: * ::: 1 2
parallel ln -s `pwd`/{1}/Chimeric.out.junction ../$DEST/{1}.Chimeric.out.junction ::: * 
parallel ln -s `pwd`/{1}/Aligned.noS.bam ../$DEST/{1}.bam ::: * 
parallel ln -s `pwd`/{1}/Aligned.noS.bam.bai ../$DEST/{1}.bam.bai ::: *
parallel ln -s `pwd`/{1}/SJ.out.tab ../$DEST/{1}.SJ.out.tab ::: *

cd ..

ls $DEST/*.bam | grep bam | grep -v bai | grep -v mate | grep -v bam_files > $DEST/bam_files.txt
ls $DEST/*.junction | grep Chimeric.out.junction | grep mate1 | grep -v fixed > $DEST/mate1
ls $DEST/*.junction | grep Chimeric.out.junction | grep mate2 | grep -v fixed > $DEST/mate2
ls $DEST/*.junction | grep Chimeric.out.junction | grep -v mate  > $DEST/samplesheet

