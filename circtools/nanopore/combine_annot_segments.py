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
    # Initialize variables
    first_name = "NONAME"
    first_col2 = None
    sum_col5 = 0
    sum_col6 = 0
    sum_col7 = 0
    sum_col8 = 0
    sum_col9 = 0

    count = 0
    internal_count = 0
    prev_col1 = None
    prev_col3 = None

    # Process each line from input
    for line in sys.stdin:
        line = line.strip()
        count += 1

        # Split line into columns
        cols = line.split("\t")
        col1, col2, col3, name, col5, col6, col7, col8, col9 = cols
        col5, col6, col7, col8, col9 = map(float, [col5, col6, col7, col8, col9])

        internal_count += 1

        if internal_count == 1:
            first_col2 = col2
            first_name = name

        if name == first_name:
            sum_col5 += col5
            sum_col6 += col6
            sum_col7 += col7
            sum_col8 += col8
            sum_col9 += col9
            prev_col1 = col1
            prev_col3 = col3

        elif count > 1:
            # Round sums to nearest integer
            round5 = int(sum_col5 + 0.5)
            round6 = int(sum_col6 + 0.5)
            round7 = int(sum_col7 + 0.5)
            round8 = int(sum_col8 + 0.5)
            round9 = int(sum_col9 + 0.5)

            # Print the accumulated result
            print(f"{prev_col1}\t{first_col2}\t{prev_col3}\t{first_name}\t{round5}\t{round6}\t{round7}\t{round8}\t{round9}")

            # Reset sums for the next group
            sum_col5 = col5
            sum_col6 = col6
            sum_col7 = col7
            sum_col8 = col8
            sum_col9 = col9
            first_col2 = col2
            first_name = name

if __name__ == "__main__":
    main()
