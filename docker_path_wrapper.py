#!/usr/bin/env python3

# Copyright (C) 2025 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import subprocess

def process_paths(paths):
    modified_paths_int = []
    for path in paths:
        if os.path.isabs(path):
            # Prefix absolute path with "/host_os"
            modified_path = os.path.join("/host_os", path.lstrip("/"))
            modified_paths_int.append(modified_path)
        # not an (absolut) path, add unmodified
        else:
            modified_paths_int.append(path)
    return modified_paths_int

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(1)

    user_paths = sys.argv[1:]
    modified_paths = process_paths(user_paths)

    modified_paths.insert(0,'circtools')

    subprocess.run(modified_paths)
