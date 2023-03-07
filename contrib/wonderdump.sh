#!/bin/bash

#
# Wonderdump is a workaround to download SRA files directly
# when fastq-dump's internet connection does not work.
# Which can happen surprisingly frequently.
#
# Usage:
#   wonderdump -X 10000 SRR1553500

# this is a slightly modified version of wonderdump
# please see http://data.biostarhandbook.com/scripts/README.html for the
# original version

set -ue

# This is where we will store the file.
SRA_DIR=./

# Make the directory if it does not exist.
mkdir -p $SRA_DIR

# All other parameters up to last.
PARAMS=${@:1:$(($# - 1))}

# The last parameter must be the SRR number.
SRR=${@:$#}

if [ -z "$SRR" ]
then
      echo "*** Please specify an SRR number"
      exit;
fi

echo "*** Getting SRR run: $SRR"

# Create the full path to the file.
SRA_FILE="$SRA_DIR/$SRR.sra"
TMP_FILE="$SRA_DIR/$SRR.tmp"

# Download only if it does not exist.
if [ ! -f $SRA_FILE ];
then
    CLEAN_SRR=${SRR:0:10}
    URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/$CLEAN_SRR/$CLEAN_SRR"
    echo "*** Downloading: $URL"
    echo "*** Saving to: $SRA_FILE"
    curl $URL > $TMP_FILE

    # Move to local file only if successful.
    mv $TMP_FILE $SRA_FILE
else
    echo "*** SRA file found: $SRA_FILE"
fi

# Run the fastq-dump.
fastq-dump $PARAMS ${SRA_FILE}
