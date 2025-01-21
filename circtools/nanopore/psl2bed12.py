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

    for line in sys.stdin:
        line = line.strip()
        count += 1

        # Skip the first 5 lines (header in PSL format)
        if count > 5:
            # Split the line into columns
            cols = line.split("\t")

            # Extract relevant columns
            col1 = cols[0]
            col9 = cols[8]
            col10 = cols[9]
            col14 = cols[13]
            col16 = int(cols[15])
            col17 = int(cols[16])
            col18 = cols[17]
            col19 = cols[18]
            col21 = cols[20]

            # Parse blockStarts and adjust them relative to col16
            block_starts = [int(x) - col16 for x in col21.split(",") if x]

            # Generate required values
            rgb = "0,0,0"
            name = f"{col10}~{col14}:{col16}-{col17}"

            # Print the output in the desired format
            print(f"{col14}\t{col16}\t{col17}\t{name}\t{col1}\t{col9}\t{col16}\t{col17}\t{rgb}\t{col18}\t{col19}\t", end="")
            print(",".join(map(str, block_starts)) + ",")

if __name__ == "__main__":
    main()