#!/usr/bin/env python3

# Copyright (C) 2024 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


def main():
    import sys

    prev_name = "NONAME"
    prev_start = None
    prev_end = None
    prev_length = None
    prev_gene = None
    prev_exon = None
    prev_est = None
    prev_intron = None

    count = 0
    internal_count = 0
    do_print = 0
    for_printer = None

    for line in sys.stdin:
        line = line.strip()
        count += 1

        if count > 1:  # Skip the header
            cols = line.split("\t")
            col1, col2, col3, col4, length, gene, exon, est, intron = cols
            length, gene, exon, est, intron = map(int, [length, gene, exon, est,
                                                        intron])
            name, region = col4.split("~")
            chrStart, end = region.split("-")
            chr_, start = chrStart.split(":")
            start, end = int(start), int(end)

            if prev_name != name:
                internal_count = 0
                if do_print == 1:
                    print(for_printer)

            internal_count += 1

            if count > 2 and internal_count == 1:
                prev_name = name
                prev_chr = chr_
                prev_start = start
                prev_end = end
                prev_length = length
                prev_gene = gene
                prev_exon = exon
                prev_est = est
                prev_intron = intron
                do_print = 0

            elif internal_count == 2:
                if prev_name == name:
                    start_min = min(prev_start, start)
                    end_max = max(prev_end, end)
                    length_sum = prev_length + length
                    gene_sum = prev_gene + gene
                    exon_sum = prev_exon + exon
                    est_sum = prev_est + est
                    intron_sum = prev_intron + intron
                    for_printer = f"{chr_}\t{start_min}\t{end_max}\t{name}\t{length_sum}\t{gene_sum}\t{exon_sum}\t{est_sum}\t{intron_sum}"
                    do_print = 1
                else:
                    do_print = 0

            elif internal_count > 2:
                do_print = 0

            if prev_name != name:
                internal_count = 1
                prev_name = name
                prev_start = start
                prev_end = end
                prev_length = length
                prev_gene = gene
                prev_exon = exon
                prev_est = est
                prev_intron = intron
                do_print = 0

    if internal_count == 2 and do_print == 1:
        print(for_printer)


if __name__ == "__main__":
    main()
