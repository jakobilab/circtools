# This module annotates circRNAs with gene names, infers library strandedness,
# and filters candidates by requiring that they originate from a single annotated gene.

import logging
import os
import re
import warnings
import HTSeq
import pysam
from .IntervalTree import IntervalTree

class CircAnnotate(object):
    def __init__(self, tmp_dir, strand=True):
        self.strand = strand
        self.tmp_dir = tmp_dir

    def clean_xenbase_name(self, name):
        """
        Extracts the standard gene symbol from a provisional Xenbase tag.
        Example: 'LOC108703370 [provisional:kcnma1]' -> 'kcnma1'
        Returns the original name if no provisional tag exists.
        """
        if not name or name == 'N/A':
            return 'N/A'

        provisional = re.search(r'\[provisional:(.*?)\]', name)
        if provisional:
            return provisional.group(1).strip()
        return name.strip()

    def parse_attributes(self, attr_data):
        """
        Safely extracts gene information from either an HTSeq attribute
        dictionary or a raw GTF attribute string.
        """
        # If HTSeq successfully parsed attributes into a dictionary
        if isinstance(attr_data, dict):
            for key in ['gene_name', 'gene', 'gene_id', 'transcript_id']:
                if key in attr_data:
                    return self.clean_xenbase_name(attr_data[key])
            return 'N/A'

        # Fallback for raw string attributes (common in some GFF3 formats)
        data = {}
        parts = [p.strip() for p in attr_data.split(';') if p.strip()]
        for p in parts:
            kv = re.match(r'^(\w+)\s+"?(.*?)"?$', p)
            if kv:
                data[kv.group(1)] = kv.group(2)

        for key in ['gene_name', 'gene', 'gene_id', 'transcript_id']:
            if key in data:
                return self.clean_xenbase_name(data[key])
        return 'N/A'

    def selectGeneGtf(self, gtf_file):
        """
        Constructs the annotation interval tree.
        Includes exons, CDS, tRNA, and rRNA to ensure proper capture of
        mitochondrial and non-coding genes.
        """
        gtf = HTSeq.GFF_Reader(gtf_file, end_included=True)
        annotation_tree = IntervalTree()
        gtf_exon = []

        valid_features = {'exon', 'CDS', 'tRNA', 'rRNA'}

        for feature in gtf:
            # Save standard exon lines for the sorted temporary file
            if feature.type == 'exon':
                gtf_exon.append(feature.get_gff_line().split('\t'))

            # Insert valid features into the interval tree
            if feature.type in valid_features:
                iv = feature.iv
                annotation_tree.insert(iv, annotation=feature)

        # Sort and write the temporary exon file
        gtf_exon_sorted = sorted(gtf_exon, key=lambda x: (x[0], int(x[3]), int(x[4])))
        gtf_exon_sorted = ['\t'.join(s) for s in gtf_exon_sorted]

        tmp_path = os.path.join(self.tmp_dir, "tmp_" + os.path.basename(gtf_file) + '.exon.sorted')
        with open(tmp_path, 'w') as new_gtf:
            new_gtf.writelines(gtf_exon_sorted)

        return annotation_tree

    def infer_library_strandedness(self, bam_file, annotation_tree, sample_size=200000):
        """
        Infers the strandedness of RNA-seq data by comparing read alignments
        against the known GTF annotation tree. Returns 'unstranded',
        'stranded_forward', 'stranded_reverse', or 'ambiguous'.
        """
        logging.info(f"Inferring library strandedness from {bam_file}...")

        sense_count = 0
        antisense_count = 0
        total_intersected = 0

        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
        except Exception as e:
            logging.warning(f"Could not open BAM for inference: {e}")
            return "unknown"

        for i, read in enumerate(bam):
            if total_intersected >= sample_size:
                break
            # Failsafe to prevent scanning a massive BAM if reads aren't mapping to genes
            if i > sample_size * 10:
                break

            if read.is_unmapped or not read.is_proper_pair:
                continue

            read_strand = '-' if read.is_reverse else '+'
            try:
                iv = HTSeq.GenomicInterval(read.reference_name, read.reference_start, read.reference_end, ".")
            except ValueError:
                continue

            out_features = []
            try:
                annotation_tree.intersect(iv, lambda x: out_features.append(x.annotation))
            except KeyError:
                continue

            if not out_features:
                continue

            gene_strand = out_features[0].iv.strand
            if gene_strand == '.':
                continue

            # ==========================================
            # --- NEW: PRINT FIRST 20 MATCHES TO CONSOLE ---
            if total_intersected < 2000:
                # Determine if it's R1 or R2
                read_type = "R1" if read.is_read1 else ("R2" if read.is_read2 else "Unk")
                # Grab the gene name for context
                gene_name = self.parse_attributes(out_features[0].attr)

                debug_msg = (f"DEBUG -> {read_type} on '{read_strand}' strand "
                             f"| Overlaps '{gene_name}' on '{gene_strand}' strand "
                             f"| ID: {read.query_name}")
                print(debug_msg)
                logging.info(debug_msg)
            # ==========================================

            total_intersected += 1

            # Determine read orientation relative to the gene
            # Determine read orientation relative to the gene
            if read.is_read1:
                if read_strand == gene_strand:
                    sense_count += 1
                else:
                    antisense_count += 1
            elif read.is_read2:
                if read_strand != gene_strand:
                    sense_count += 1
                else:
                    antisense_count += 1

        bam.close()

        if total_intersected < 100:
            logging.warning("Not enough mapped reads overlapping genes to infer strandedness.")
            return "unknown"

        fraction_sense = sense_count / total_intersected
        fraction_antisense = antisense_count / total_intersected

        logging.info(f"Strand inference based on {total_intersected} reads:")
        logging.info(f"Sense: {fraction_sense:.2f} | Antisense: {fraction_antisense:.2f}")

        if fraction_sense > 0.8:
            return "stranded_forward"
        elif fraction_antisense > 0.8:
            return "stranded_reverse"
        elif 0.4 <= fraction_sense <= 0.6:
            return "unstranded"
        else:
            return "ambiguous"

    def annotate_one_interval(self, interval, annotation_tree, what='gene'):
        """
        Intersects a single genomic interval with the annotation tree.
        Extracts either the gene name or the region type.
        """
        # Force strand to '.' to allow mapping if the library is unstranded
        if not self.strand:
            interval.strand = '.'

        out_features = []
        annotation_tree.intersect(interval, lambda x: out_features.append(x.annotation))

        collect = set()
        for feat in out_features:
            if what == 'gene':
                name = self.parse_attributes(feat.attr)
                collect.add(name)
            else:
                collect.add(feat.type)

        if 'N/A' in collect and len(collect) > 1:
            collect.remove('N/A')

        genes = ','.join(sorted(collect))
        return genes if genes else 'not_annotated'

    def annotate(self, circfile, annotation_tree, output):
        """
        Iterates through the circRNA BED file and annotates each candidate
        with its corresponding host gene.
        """
        with open(circfile, 'r') as tmpcirc:
            line = tmpcirc.readline()
            if line and len(line.split('\t')) != 6:
                warnings.warn('Input circRNA file is not in the desired bed6 format!')

        out = open(output, 'w')
        circ_regions = HTSeq.BED_Reader(circfile)
        for circ in circ_regions:
            annotation = self.annotate_one_interval(circ.iv, annotation_tree, what='gene')
            out.write('\t'.join([
                circ.iv.chrom,
                str(circ.iv.start),
                str(circ.iv.end),
                annotation,
                str(int(circ.score)),
                circ.iv.strand
            ]) + '\n')
        out.close()

    def auto_annotate(self, circfile, gtf_file, bam_file, output):
        """
        Convenience wrapper: builds the interval tree, checks the BAM for
        strandedness, adjusts parameters, and annotates the circRNAs.
        """
        logging.info("Building GTF annotation tree...")
        tree = self.selectGeneGtf(gtf_file)

        inferred_type = self.infer_library_strandedness(bam_file, tree)

        if inferred_type == "unstranded":
            logging.info("Unstranded library detected. Disabling strict strand overlap requirement.")
            self.strand = False
        elif inferred_type in ["stranded_forward", "stranded_reverse"]:
            logging.info(f"Directional library ({inferred_type}) detected. Enforcing strict strand overlap.")
            self.strand = True
        else:
            logging.warning("Strandedness ambiguous. Defaulting to strict strand enforcement.")

        logging.info(f"Annotating {circfile}...")
        self.annotate(circfile, tree, output)

    def uniqstring(self, strings, sep=','):
        """
        Helper function to remove duplicate annotations and sort them.
        """
        if not strings or strings == 'not_annotated':
            return 'not_annotated'
        parts = set(strings.split(sep))
        if '.' in parts: parts.remove('.')
        return sep.join(sorted(parts))

    def annotateregions(self, circfile, annotation_tree, output):
        """
        Annotates the left and right boundaries of the circRNA to determine
        if the junctions fall within exons, introns, or intergenic regions.
        """
        circ = open(circfile, 'r').readlines()
        with open(output, 'w') as new_CircCoordinates:
            new_CircCoordinates.write('Chr\tStart\tEnd\tGene\tJunctionType\tStrand\tStart-End Region\tOverallRegion\n')
            for line in circ:
                line_split = line.split('\t')

                # Create 1-bp intervals to test left and right splice boundaries
                iv_left = HTSeq.GenomicInterval(line_split[0], int(line_split[1]), int(line_split[1]) + 1, line_split[5].strip())
                iv_right = HTSeq.GenomicInterval(line_split[0], int(line_split[2]) - 1, int(line_split[2]), line_split[5].strip())

                left_ann = self.annotate_one_interval(iv_left, annotation_tree, what='region')
                right_ann = self.annotate_one_interval(iv_right, annotation_tree, what='region')

                overall = self.uniqstring(left_ann + ',' + right_ann)

                left_label = self.readRegionAnnotate(left_ann)
                right_label = self.readRegionAnnotate(right_ann)

                new_CircCoordinates.write(line.strip() + '\t' + left_label + '-' + right_label + '\t' + overall + '\n')

    def readRegionAnnotate(self, annotations):
        """
        Classifies the genomic region based on the intersected features.
        """
        if 'exon' in annotations or 'CDS' in annotations:
            return 'exon'
        elif 'not_annotated' in annotations:
            return 'intergenic'
        else:
            return 'intron'

    def filtbygene(self, circ2filter, output):
        """
        Filters the circRNA candidates, keeping only those that originate
        from exactly one annotated gene. Discards intergenic and multi-gene hits.
        """
        out = open(output, 'w')
        with open(circ2filter, 'r') as circ:
            for line in circ:
                tmp = line.split('\t')
                genes = tmp[3].split(',')
                # Only keep candidates mapped to a single, valid gene
                if len(genes) == 1:
                    out.write(line)
        out.close()