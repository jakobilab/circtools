#!/usr/bin/env python2

import collections
import logging
import re


class Findcirc(object):
    # Initialize some parameters
    def __init__(self, endTol, minL, maxL, library_type="stranded_reverse"):
        # self.strand = strand
        self.endTol = endTol
        # self.output = output
        self.maxL = int(maxL)
        self.minL = int(minL)
        self.library_type = library_type

    # self.findcirc(Chim_junc)

    def cigarGenomicDist(self, cig):
        C = re.findall('[a-zA-Z]', cig)
        L = re.findall('\-?[0-9]+', cig)
        n = len(L)
        g = 0
        for i in range(0, n):
            if C[i] != 'S' and C[i] != "I":
                g = g + int(L[i])
        return g

    def printcircline(self, Chim_junc, output):
        junctionfile = open(Chim_junc, 'r')
        outfile = open(output, 'w')
        for line in junctionfile:
            L = line.split('\t')
            if L[0] == "chr_donorA":
                continue
            if int(L[6]) >= 0 and L[0] == L[3] and L[2] == L[5] and (
                    (L[2] == '-' and int(L[4]) > int(L[1]) and self.minL < (int(L[4]) - int(L[1])) < self.maxL) or (
                    L[2] == '+' and int(L[1]) > int(L[4]) and self.minL < (
                    int(L[1]) - int(L[4])) < self.maxL)):
                if (L[2] == '+' and (int(L[10]) + self.endTol) > int(L[4]) and (
                        int(L[12]) + self.cigarGenomicDist(L[13]) - self.endTol) <= int(L[1])) or (
                        L[2] == '-' and (int(L[12]) + self.endTol) > int(L[1]) and (
                        int(L[10]) + self.cigarGenomicDist(L[11]) - self.endTol) <= int(L[4])):
                    outfile.write(line)
        outfile.close()
        junctionfile.close()

    def sepDuplicates(self, Chim_junc, duplicates, nonduplicates):
        junctionfile = open(Chim_junc, 'r')
        dup = open(duplicates, 'w')
        nondup = open(nonduplicates, 'w')

        reads = []
        lines = []  # A list of all lines
        suffice = False
        for line in junctionfile:
            line_split = line.split('\t')
            if not suffice:
                if len(line_split[9].split('.')[-1]) == 1:
                    suffice = True
            if suffice:
                readname = '.'.join(
                    (line_split[1], line_split[2], line_split[4], '.'.join(line_split[9].split('.')[:-1])))
            else:
                readname = '.'.join((line_split[1], line_split[2], line_split[4], line_split[9]))
            reads.append(readname)
            lines.append(line)

        for indx, read in enumerate(reads):
            if reads.count(read) == 2:
                dup.write(lines[indx])
            elif reads.count(read) > 2:
                print('Read %s has more than 2 count.' % read)
                try:
                    logging.warning('Read %s has more than 2 count.' % read)
                except NameError:
                    pass
            else:
                nondup.write(lines[indx])

        junctionfile.close()
        dup.close()
        nondup.close()

    def smallcirc(self, duplicates, output, strand=True):
        dup = open(duplicates).readlines()
        outfile = open(output, 'w')

        collect = []
        for line in dup:
            L = line.split('\t')
            if L[2] == '+':
                identifier = (L[0], L[1], L[3], L[4], L[9])
            elif L[2] == '-':
                identifier = (L[0], L[4], L[3], L[1], L[9])
            if identifier in collect:
                # Dynamic strand selection based on library type
                if self.library_type == "stranded_forward":
                    target_strand = L[2]
                else:
                    target_strand = '-' if L[2] == '+' else '+'

                if L[2] == '+':
                    if strand:
                        score = '0' if L[6] == '0' else str(3 - int(L[6]))
                        res = [L[0], str(int(L[4]) + 1), str(int(L[1]) - 1), target_strand, score, L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
                    else:
                        res = [L[0], str(int(L[4]) + 1), str(int(L[1]) - 1), L[2], L[6], L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
                if L[2] == '-':
                    if strand:
                        score = '0' if L[6] == '0' else str(3 - int(L[6]))
                        res = [L[0], str(int(L[1]) + 1), str(int(L[4]) - 1), target_strand, score, L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
                    else:
                        res = [L[0], str(int(L[1]) + 1), str(int(L[4]) - 1), L[2], L[6], L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
            else:
                collect.append(identifier)

        outfile.close()

    def findcirc(self, Chim_junc, output, strand=True):
        junctionfile = open(Chim_junc, 'r')
        outfile = open(output, 'w')
        linecnt = 1
        for line in junctionfile:
            L = line.split('\t')
            linecnt = linecnt + 1

            if len(L) < 14:
                print(("WARNING: File " + str(Chim_junc) + ", line " + str(linecnt) + " does not contain all features."))
                print(("WARNING: " + str(Chim_junc) + " is probably corrupt."))
            if L[0] == "chr_donorA":
                continue
            if int(L[6]) >= 0 and L[0] == L[3] and L[2] == L[5] and (
                    (L[2] == '-' and int(L[4]) > int(L[1]) and self.minL < (int(L[4]) - int(L[1])) < self.maxL) or (
                    L[2] == '+' and int(L[1]) > int(L[4]) and self.minL < (
                    int(L[1]) - int(L[4])) < self.maxL)):

                # Determine output strand based on library geometry
                if self.library_type == "stranded_forward":
                    pos_match_strand = '+'
                    neg_match_strand = '-'
                else: # Default behavior (dUTP / stranded_reverse)
                    pos_match_strand = '-'
                    neg_match_strand = '+'

                if (L[2] == '+' and (int(L[10]) + self.endTol) > int(L[4]) and (
                        int(L[12]) + self.cigarGenomicDist(L[13]) - self.endTol) <= int(L[1])):
                    if strand:
                        score = '0' if L[6] == '0' else str(3 - int(L[6]))
                        res = [L[0], str(int(L[4]) + 1), str(int(L[1]) - 1), pos_match_strand, score, L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
                    else:
                        res = [L[0], str(int(L[4]) + 1), str(int(L[1]) - 1), L[2], L[6], L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
                if L[2] == '-' and (int(L[12]) + self.endTol) > int(L[1]) and (
                        int(L[10]) + self.cigarGenomicDist(L[11]) - self.endTol) <= int(L[4]):
                    if strand:
                        score = '0' if L[6] == '0' else str(3 - int(L[6]))
                        res = [L[0], str(int(L[1]) + 1), str(int(L[4]) - 1), neg_match_strand, score, L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
                    else:
                        res = [L[0], str(int(L[1]) + 1), str(int(L[4]) - 1), L[2], L[6], L[7], L[8]]
                        outfile.write(('\t').join(res) + '\n')
        outfile.close()
        junctionfile.close()


class Sort(object):
    def __init__(self):
        pass

    def count(self, sortedlist, strand=True):
        cnt = collections.Counter()
        tmp_count = []
        for itm in sortedlist:
            if strand:
                circs = (itm[0], itm[1], itm[2], itm[3])
            elif not strand:
                circs = (itm[0], itm[1], itm[2])
            else:
                print("Please specify correct strand information.")
            cnt[circs] += 1
            itm.append(str(cnt[circs]))
            tmp_count.append([itm[0], itm[1], itm[2], '.', itm[7], itm[3], itm[4], itm[5], itm[6]])
        tmp_count.reverse()
        tmp_count_dedup = []
        lines_seen = set()
        for itm in tmp_count:
            if strand:
                tmp_itm = tuple((itm[0], itm[1], itm[2], itm[5]))
            elif not strand:
                tmp_itm = tuple((itm[0], itm[1], itm[2]))
            if tmp_itm not in lines_seen:
                tmp_count_dedup.append(itm)
                lines_seen.add(tmp_itm)
        tmp_count_dedup.reverse()
        return tmp_count_dedup

    def sort_count(self, findcircOut, output, strand=True):
        file_list = open(findcircOut, 'r').readlines()
        output = open(output, 'w')
        tmp_sort = []
        for itm in file_list:
            line_tmp = itm.strip().split('\t')
            tmp_sort.append(line_tmp)
        tmp_sorted = sorted(tmp_sort, key=lambda x: (x[0], int(x[1]), int(x[2]), x[5]))
        sorted_count = self.count(tmp_sorted, strand=strand)
        output.writelines('\t'.join(j) + '\n' for j in sorted_count)
        output.close()