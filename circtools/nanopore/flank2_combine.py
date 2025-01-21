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
    first_seq = ""

    for line in sys.stdin:
        line = line.strip()
        count += 1

        col1, col2 = line.split("\t")
        name, chr_, start, end, strand = col1.split("~")
        start, end = int(start), int(end)
        length = end - start

        if count == 1:
            first_seq = col2
        elif count == 2:
            count = 0
            splice_site = f"{first_seq}{col2}"
            if splice_site == "AGGT":
                genome_strand = "+"
            elif splice_site == "ACCT":
                genome_strand = "-"
            else:
                genome_strand = "Unknown"

            print(f"{chr_}\t{start}\t{end}\t{name}\t{length}\t{genome_strand}\t{splice_site}")

if __name__ == "__main__":
    main()