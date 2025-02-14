#!/bin/bash

# This file contains code that is derived from MIT-licensed software:
#   Original code by Morten T Venø <morten.veno@omiics.com>,
#   licensed under the MIT License.
#
# The MIT License:
# ---------------------------------------------------------------------
# See LICENSE_MIT for the full MIT license text.
# ---------------------------------------------------------------------
# Modifications and additions to this file are licensed under:
# The GNU General Public License (GPL), version 3 or later.
#
# Copyright (C) 2025 Tobias Jakobi
#
# See the LICENSE file in the root directory for the full terms of the GPL.
#


# This is version 5.5 of the Nanopore circRNA datection pipeline
# This version uses blat to map reads.  NanoFilt is used to trim low quality reads (q 7). Parallel BLAT is used to map nanopore reads to genome.

if [ $# -ne 8 ]; then
    echo "This script is called directly by the circtools nanopore pipeline and not meant for direct user interaction."
    exit 1
fi

data_folder=$1
sample=$2
genome=$3
reference_path=$4
scriptFolder=$5
output_path=$6
threads=$7
sample_ext=$8

mkdir "$output_path"
cd "$output_path" || exit

fa=$reference_path/$genome/genome.fa
mRNA=$reference_path/$genome/refFlat.csv.unique.bed
exon=$reference_path/$genome/refFlat.csv.merged.bed
single_exon=$reference_path/$genome/refFlat.csv.sort.bed
est=$reference_path/$genome/est.bed
genomeSize=$reference_path/$genome/genome.chrom.sizes
circRNA_prefix="$genome"_circ_


# Temp folder
temp_sort=$(mktemp -d /tmp/foo.XXXXXXXXX)

date

echo "Sample: "$sample
echo
echo "Number of raw reads before NanoFilt -q 7 -l 250"
zcat $data_folder/$sample.$sample_ext | wc -l | awk '{print $1/4}'

echo
date
echo "NanoFilt to remove reads under quality 7 and conversion to fasta"
zcat $data_folder/$sample.$sample_ext | NanoFilt -q 7 -l 250 | sed -n '1~4s/^@/>/p;2~4p' > $sample.fa
echo
date
echo "Number of filtered reads after NanoFilt -q 7 -l 250"
cat $sample.fa | wc -l | awk '{print $1/2}'
echo
date
echo "Mapping with pblat - parallelized blat with multi-threads support (http://icebert.github.io/pblat/)"
echo "lower case sequences in the genome file are masked out"
echo "Showing a dot for every 50k sequences processed"
pblat  -threads=$threads -trimT -dots=50000 -mask=lower $fa $sample.fa $sample.psl
echo "Blat done"
date

# Remove non-standard chromosomes
cat $sample.psl | grep -v "_random" | grep -v "_hap" | grep -v "chrUn_" > $sample.temp.psl
rm $sample.psl
mv $sample.temp.psl $sample.psl
cat $sample.psl | awk '{print $10}' | sort | uniq -c | sort -nrk 1,1 > mappings_per_read.txt
cat mappings_per_read.txt | awk '{print $1}' | uniq -c | sort -nrk 1,1 > $sample.histogram_number_of_genomic_hits_per_read.txt

cat $sample.psl | python3 $scriptFolder/psl2bed12.py |  bedtools sort > $sample.psl.bed

echo
date
echo "Getting group numbers"
cat $sample.psl | python3 $scriptFolder/blat_output_processing_v3.py > $sample.scan.psl

# Read fragment numbers
echo
date
echo "The different groups, numbers of read fragments:"
cat $sample.scan.psl | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6
cat $sample.scan.psl | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6 > $sample.scan.groupNumbers.fragments.txt
# Read numbers
echo
echo "The different groups, numbers of unique reads:"
cat $sample.scan.psl | awk '{print $10,$NF}' | sort | uniq | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6
cat $sample.scan.psl | awk '{print $10,$NF}' | sort | uniq | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 | head -6 > $sample.scan.groupNumbers.reads.txt


# Get circRNA reads
head -5 $sample.scan.psl > $sample.scan.circRNA.psl
grep circRNA $sample.scan.psl >> $sample.scan.circRNA.psl
## converting psl to bed12
cat $sample.scan.circRNA.psl | python3 $scriptFolder/psl2bed12.py |  bedtools sort > $sample.scan.circRNA.psl.bed



#echo
#echo "Making bam and bigWig (.bw) files for use in genome browsers"
## All Blat mapped reads
#bedtools bedtobam -i $sample.psl.bed -bed12 -g $genomeSize > $sample.bam
#samtools sort $sample.bam > $sample.sort.bam
#samtools index $sample.sort.bam
#rm $sample.bam
#bamCoverage --binSize 1 --numberOfProcessors $threads -b $sample.sort.bam -o $sample.bw
#
## Blat mapped BSJ spanning reads
#bedtools bedtobam -i $sample.scan.circRNA.psl.bed -bed12 -g $genomeSize > $sample.circRNA.bam
#samtools sort $sample.circRNA.bam > $sample.circRNA.sort.bam
#samtools index $sample.circRNA.sort.bam
#rm $sample.circRNA.bam
#bamCoverage --binSize 1 --numberOfProcessors $threads -b $sample.circRNA.sort.bam -o $sample.circRNA.bw



echo
date
echo "outputting Potential_multi-round_circRNA"
# Get Potential_multi-round_circRNA
head -5 $sample.scan.psl > $sample.scan.Potential_multi-round_circRNA.psl
grep Potential_multi-round_circRNA $sample.scan.psl >> $sample.scan.Potential_multi-round_circRNA.psl
## converting psl to bed12
cat $sample.scan.Potential_multi-round_circRNA.psl | python3 $scriptFolder/psl2bed12.py |  bedtools sort > $sample.scan.Potential_multi-round_circRNA.psl.bed
bedtools bedtobam -i $sample.scan.Potential_multi-round_circRNA.psl.bed -bed12 -g $genomeSize > $sample.scan.Potential_multi-round_circRNA.bam
samtools sort $sample.scan.Potential_multi-round_circRNA.bam > $sample.scan.Potential_multi-round_circRNA.sort.bam
samtools index $sample.scan.Potential_multi-round_circRNA.sort.bam
# Making a special file to show how many rounds each read takes
cat $sample.scan.Potential_multi-round_circRNA.psl.bed | sed 's/~/\t/g' | awk 'OFS="\t"{print $4,$2,$3,$1,$6,$7}' |  bedtools sort | awk 'OFS="\t"{print $1"~"$4,$2,$3,$4,$5,$6}' | bedtools merge > $sample.scan.Potential_multi-round_circRNA.psl.merge.bed
cat $sample.scan.Potential_multi-round_circRNA.psl.bed | sed 's/~/\t/g' | awk 'OFS="\t"{print $4,$2,$3,$1,$6,$7}' |  bedtools sort | awk 'OFS="\t"{print $1"~"$4,$2,$3,$4,$5,$6}' | bedtools coverage -counts -a $sample.scan.Potential_multi-round_circRNA.psl.merge.bed -b - > temp.$sample.multi-round.count.txt
printf "#Chr\tStart\tEnd\tRead_name\tNumber_of_rounds\tOverlapping_gene\n" > $sample.scan.Potential_multi-round_circRNA.psl.annot.bed
cat temp.$sample.multi-round.count.txt | sed 's/~/\t/g' | awk 'OFS="\t"{print $2,$3,$4,$1,$5}' |  bedtools sort | bedtools map -c 4 -o distinct -a - -b $mRNA >> $sample.scan.Potential_multi-round_circRNA.psl.annot.bed
cat $sample.scan.Potential_multi-round_circRNA.psl.annot.bed | awk '{print $NF}' | sort | uniq -c | sort -nrk 1,1 > $sample.scan.Potential_multi-round_circRNA.psl.annot.count.txt
cat $sample.scan.Potential_multi-round_circRNA.psl.annot.bed | awk '{print $4}' | grep -v Read_name > $sample.temp.read_names
#grep --no-group-separator -A1 -f $sample.temp.read_names $sample.fa > $sample.Potential_multi-round_circRNA.fa
samtools faidx $sample.fa
echo "Using samtools to make Potential_multi-round_circRNA fasta files for sample: "$sample.fa
xargs samtools faidx $sample.fa < $sample.temp.read_names > $sample.Potential_multi-round_circRNA.fa
rm temp.$sample.multi-round.count.txt $sample.temp.read_names

echo
date
echo "Annotating circRNAs and converting to bed6 format"

### Adding split
echo "  First splitting bed file in 8 to optimize run time"
bedtools split -n 8 -p $sample.scan.circRNA.psl -a simple -i $sample.scan.circRNA.psl.bed

## Very time consuming step for large datasets... Converting to bed6 (turning each block from bed12 into a single bed entry). Header is removed using awk
cat $sample.scan.circRNA.psl.00001.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00001.bed &
cat $sample.scan.circRNA.psl.00002.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00002.bed &
cat $sample.scan.circRNA.psl.00003.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00003.bed &
cat $sample.scan.circRNA.psl.00004.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00004.bed &
cat $sample.scan.circRNA.psl.00005.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00005.bed &
cat $sample.scan.circRNA.psl.00006.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00006.bed &
cat $sample.scan.circRNA.psl.00007.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00007.bed &
cat $sample.scan.circRNA.psl.00008.bed | bedtools bed12tobed6 -i stdin | bedtools annotate -names exons est genes -both -i stdin -files $mRNA $exon $est | awk 'NR>1 {print $0}'> $sample.scan.circRNA.psl.annot.00008.bed &

wait

cat $sample.scan.circRNA.psl.annot.00001.bed $sample.scan.circRNA.psl.annot.00002.bed $sample.scan.circRNA.psl.annot.00003.bed $sample.scan.circRNA.psl.annot.00004.bed $sample.scan.circRNA.psl.annot.00005.bed $sample.scan.circRNA.psl.annot.00006.bed $sample.scan.circRNA.psl.annot.00007.bed $sample.scan.circRNA.psl.annot.00008.bed > $sample.scan.circRNA.psl.annot.bed
rm $sample.scan.circRNA.psl.annot.0*.bed

### Done with split



printf "#chr\tstart\tend\tread_name\tread_length\tgene_coverage\texon_coverage\tEST_coverage\tintron_coverage\n" > $sample.scan.circRNA.psl.annot.txt
cat $sample.scan.circRNA.psl.annot.bed | sort -k 4,4 -T $temp_sort | awk 'OFS="\t"{print $1,$2,$3,$4,$3-$2,($3-$2)*$8,($3-$2)*$10,($3-$2)*$12,($3-$2)*$8-($3-$2)*$10}' | python3 $scriptFolder/combine_annot_segments.py >> $sample.scan.circRNA.psl.annot.txt
cat $sample.scan.circRNA.psl.annot.txt | python3 $scriptFolder/make_circRNAs_from_annot.txt.py > $sample.scan.circRNA.psl.annot.combine.txt

echo
date
echo "Refining circRNA edges based annotated exon boundaries and annotated circRNAs"
# Making a unique list of circRNAs
cat $sample.scan.circRNA.psl.annot.combine.txt | awk 'NR>1,OFS="\t"{print $1,$2,$2+1,$4"~"$5"~"$6"~"$7"~"$8"~"$9}' | bedtools sort | uniq > temp_start
cat $sample.scan.circRNA.psl.annot.combine.txt | awk 'NR>1,OFS="\t"{print $1,$3-1,$3,$4"~"$5"~"$6"~"$7"~"$8"~"$9}' | bedtools sort | uniq > temp_end
# bedtools sort -i $single_exon > exon_ref
# Prints the start and end position of closest exon. In special cases where the circRNA is produced far inside an annoteted exon, such as occurs for Malat1, are filtered away.
bedtools closest -t first -d -header -a temp_start -b $single_exon -nonamecheck | awk 'OFS="\t"{if ($2 - $6 < 31) print $1,$6,$7,$4,$10,$11}' | awk 'OFS="\t"{if($2 > -1) print $0 }' > temp_start.exon
bedtools closest -t first -d -header -a temp_end -b $single_exon -nonamecheck | awk 'OFS="\t"{if ($7 - $3 < 31) print $1,$6,$7,$4,$10,$11}' | awk 'OFS="\t"{if($2 > -1) print $0 }' > temp_end.exon
cat temp_start.exon temp_end.exon | sort -k 4,4 > temp_edge_exon
bedtools groupby -g 4 -c 1,2,3,5,6,4 -o distinct,min,max,distinct,max,count -i temp_edge_exon | awk 'OFS="\t"{print $2,$3,$4,$1,$6,$5,$7}' > temp_exon-ends0
# Allowing only edges that are formed from 2 read segments:
cat temp_exon-ends0 | awk 'OFS="\t"{if ($7 == 2) print $1,$2,$3,$4,$5,$6}' > temp_exon-ends
## Exon match: Max 30 bp distance, correct strand
grep -v "+,-" temp_exon-ends | awk 'OFS="\t"{if($5 < 31) print $0 }' > $sample.scan.circRNA.psl.annot.combine.correct.bed
cat $sample.scan.circRNA.psl.annot.combine.correct.bed | awk 'OFS="\t"{print $1,$2,$3,"exon_match",0,$6}' | sort -nk 3,3 | sort -nk 2,2 | sort -k 1,1 | uniq  |  bedtools sort > base_list_exon-match.bed
cat $sample.scan.circRNA.psl.annot.combine.correct.bed | awk 'OFS="\t"{print $1,$2,$3,$5,$6,$4}' | sed 's/~/\t/g' | awk 'OFS="\t"{print $1,$2,$3,$6,$4,$5,$7,$8,$9,$10,$11}' |  bedtools sort > $sample.scan.circRNA.psl.annot.combine.correct.full.bed
## No exon match: Over 30 bp distance or segments on different strands
grep "+,-" temp_exon-ends | awk 'OFS="\t"{ print $4 }' > temp_exon-ends_nohit
cat temp_exon-ends | awk 'OFS="\t"{if($5 > 30) print $4 }' >> temp_exon-ends_nohit
cat temp_exon-ends_nohit | sed 's/~/\t/g' | awk 'OFS="\t"{print $1}' | sort | uniq > temp_exon-ends_nohit_uniq

#rm temp_start temp_end exon_ref temp_start.exon temp_end.exon temp_edge_exon temp_exon-ends temp_exon-ends_nohit temp_exon-ends0

### For the base_list_exon-match.bed file
# finding host gene: Any overlap of the circRNA with an annotated refSeq gene on the SENSE strand
bedtools map -s -c 4 -o distinct -a base_list_exon-match.bed -b $mRNA -nonamecheck > base_list_exon-match.temp

awk '{print $0"\t.\t.\t."}' base_list_exon-match.temp > base_list_exon-match.temp4.bed

# print out 3 dots for each db, do we do not loose anything

# Mapping stuff on to the base list
bedtools map -f 1.0 -F 1.0 -c 4,7,8,9,10,11,5,5,5 -o count,mean,mean,mean,mean,mean,min,max,mean -a base_list_exon-match.temp4.bed -b $sample.scan.circRNA.psl.annot.combine.correct.full.bed -nonamecheck | awk 'OFS="\t"{print $1,$2,$3,$4,$11,$6,$7,$8,$9,$10,$12,$13,$14,$15,$16,$17,$18,$19}' | sort -nrk 5,5 > base_list_exon-match.annot.prefilter.bed

echo
cat base_list_exon-match.annot.prefilter.bed | grep -v chrM | grep -v Rn45s > $sample.base_list_exon-match.annot.bed

rm base_list_exon-match.temp base_list_exon-match.temp4.bed base_list_exon-match.annot.prefilter.bed

# For v 5.5 I increased the stringency. See below
### for the reads that do not match exons I check for similarity to circBase, circAtlas or CIRCpedia circRNA. If this is found, the annotated circRNA entry defines boundaries.
cat $sample.scan.circRNA.psl.annot.combine.txt |  bedtools sort > $sample.scan.circRNA.psl.annot.combine.sort.txt
# finding host gene: Any overlap of the circRNA with an annotated refSeq gene on either
bedtools intersect -wao -a $sample.scan.circRNA.psl.annot.combine.sort.txt -b $mRNA -nonamecheck | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $NF,$0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{$NF=""; print $0}' | awk 'OFS="\t"{print $0,$1}' | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' | sed s'/\t\t/\t/g' > $sample.scan.circRNA.psl.annot.combine.sort.temp

awk '{print $0"\t.\t.\t."}' $sample.scan.circRNA.psl.annot.combine.sort.temp > $sample.scan.circRNA.psl.annot.combine.circID.bed

# Getting only the no_match reads:
grep -Fwf temp_exon-ends_nohit_uniq $sample.scan.circRNA.psl.annot.combine.circID.bed |  bedtools sort > no_exon_match_reads.bed


echo
date
echo "Generating internal_circRNA_name and outputting candidate circRNA list"
# Combine the positive hits
cat $sample.base_list_exon-match.annot.bed  | sort -nrk 5,5 > temp.circ.hits

printf "chr\tstart\tend\tdescription\tBSJ_reads\tstrand\tgene\tcircBase_ID\tcircAtlas_ID\tCIRCpedia_ID\tmean_read_coverage\tmean_gene_coverage\tmean_exon_coverage\tmean_EST_coverage\tmean_intron_coverage\tmin_exon_adjust\tmax_exon_adjust\tmean_exon_adjust\n" > $sample.circRNA_candidates.annotated.bed
cat temp.circ.hits >> $sample.circRNA_candidates.annotated.bed

        count=0
        while read -r line
        do
                       if [ $count == 0 ]; then
                               echo "internal_circRNA_name" > circRNA_name.temp
                       else
                               name=$(echo $line | awk '{print $7}')
                               #antisense_name=$(echo $line | awk '{print $8}')
                               if [ $name = "." ] || [ $name = "NA" ]; then
                               #        if [ $antisense_name = "." ] || [ $antisense_name = "NA" ]; then
                                               circ_host="intergenic"
                               #        else
                               #                circ_host=$antisense_name"-AS"
                               #        fi
                               else
                                       circ_host=$name
                               fi
                               echo `printf $circRNA_prefix%04d $count`"_$circ_host" >> circRNA_name.temp
                       fi
               count=$(expr $count + 1)
        done < $sample.circRNA_candidates.annotated.bed
	paste circRNA_name.temp $sample.circRNA_candidates.annotated.bed > $sample.circRNA_candidates.annotated.txt


#### for the reads that do not match exons and also does not have 99% similarity to known circRNAs
cat no_exon_match_reads.bed | uniq | grep "\.[[:space:]]\.[[:space:]]\." |  bedtools sort > no_exon_no_circRNA.bed


## Delete temp files
rm temp.circ.hits $sample.scan.circRNA.psl.annot.combine.sort.temp temp_exon-ends_nohit_uniq $sample.scan.circRNA.psl.annot.combine.correct.full.bed circRNA_name.temp $sample.scan.circRNA.psl.annot.combine.correct.bed
rm -r $temp_sort
rm temp* $sample.scan.Potential_multi-round_circRNA.bam $sample.scan.Potential_multi-round_circRNA.psl.bed $sample.scan.Potential_multi-round_circRNA.psl
#rm $sample.scan.circRNA.psl.annot.combine.circID.bed $sample.scan.circRNA.psl.annot.combine.sort.txt $sample.base_list_exon-match.annot.bed
rm $sample.scan.Potential_multi-round_circRNA.psl.merge.bed $sample.scan.Potential_multi-round_circRNA.sort.bam.bai
rm $sample.scan.Potential_multi-round_circRNA.sort.bam
#rm $sample.scan.circRNA.bam $sample.scan.circRNA.psl.annot.bed $sample.scan.circRNA.psl.annot.txt $sample.scan.circRNA.psl $sample.sort.bam $sample.scan.circRNA.psl.bed
#rm $sample.sort.bam.bai $sample.psl.bed $sample.fa.fai $sample.psl $sample.fa no_exon_no_circRNA.bed
rm base_list_exon-match.bed no_exon_match_reads.bed
#rm base_list_no-exon.cirBaseID.annot.prefilter.bed
rm mappings_per_read.txt $sample.scan.psl #$sample.base_list_no-exon.cirBaseID.annot.bed

echo
date
echo
echo
echo Done

