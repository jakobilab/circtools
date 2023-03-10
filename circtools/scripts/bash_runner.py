#!/usr/bin/env python3

# Copyright (C) 2023 Tobias Jakobi
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

# runner idea from https://stackoverflow.com/a/73649412/1900920

import subprocess
import circtools as _

base_path = _.__path__[0]


def _run(bash_script):
    return subprocess.call(bash_script, shell=True)


def wonderdump():
    return _run("/" + base_path + "/scripts/wonderdump")

