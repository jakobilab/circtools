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

import sys


def main():
    count = 0
    internal_count = 0
    prev_type = None
    first_line = ""
    first_name = "NONAME"  # Default name

    for line in sys.stdin:
        count += 1

        # First 5 lines in PSL format are headers, print them directly
        if count <= 5:
            print(line, end="")
            continue

        line = line.strip()

        # Parse the current line
        fields = line.split("\t")
        (col1, col2, col3, col4, col5, col6, col7, col8, strand, name,
         col11, read_start, read_end, chr_, col15, start, end, *rest) = fields

        start, end = int(start), int(end)
        read_start, read_end = int(read_start), int(read_end)

        internal_count += 1

        # Process fragments from the same read
        if name == first_name:
            min_start = min(first_start, start)
            max_end = max(first_end, end)
            abs_distance = max_end - min_start

            # Has to be on the same chromosome and within a megabase
            if first_chr == chr_ and abs_distance < 1000000:
                if internal_count > 1:
                    if first_strand != strand:
                        type_ = "Not_same_strand"
                    elif strand == "+":
                        if first_read_end > read_start and (read_end - first_read_start) < 50:

                            if first_start < start:
                                type_ = "circRNA"
                            elif first_start > start:
                                type_ =  "linear"

                        elif read_end > first_read_start and (first_read_end - read_start) < 50:

                            if first_start > start:
                                type_ = "circRNA"
                            elif first_start < start:
                                type_ =  "linear"
                        else:
                            type_ = "ambiguous"
                    elif strand == "-":
                        end_to_start = abs(first_read_end - read_start)
                        start_to_end = abs(first_read_start - read_end)

                        if first_read_end > read_start and (read_end - first_read_start) < 50:
                            if first_start > start:
                                type_ = "circRNA"
                            elif first_start < start:
                                type_ =  "linear"

                        elif read_end > first_read_start and (first_read_end - read_start) < 50:

                            if first_start < start:
                                type_ = "circRNA"
                            elif first_start > start:
                                type_ =  "linear"

                        else:
                            type_ = "ambiguous"

                    # Print outputs based on fragment counts
                    if internal_count == 2:
                        print(f"{first_line}\t{type_}")
                        print(f"{line}\t{type_}")
                        prev_type = type_
                    elif internal_count > 2 and prev_type != "circRNA":
                        print(f"{line}\t{type_}")
        else:
            if count > 6:
                if internal_count == 2:
                    type_ = "linear_1_fragment"
                    print(f"{first_line}\t{type_}")

            # Start reading a new fragment
            first_line = line
            first_strand = strand
            first_name = name
            first_read_start, first_read_end = read_start, read_end
            first_chr = chr_
            first_start, first_end = start, end
            internal_count = 1
            type_ = "Potential_multi-round_circRNA"


if __name__ == "__main__":
    main()
