#!/bin/bash

# Copyright (C) 2017 Tobias Jakobi
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

function install_bedtools {

    if [[ "$OSTYPE" == "linux-gnu" ]]; then
        wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary -O /usr/bin/bedtools
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        curl -O https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary -o /usr/bin/bedtools
    else
        echo "Sorry, this OS type not supported. Please contact circtools@dieterichlab.org for help."
    fi
}

BEDTOOLS=`which bedtools`

if [ $BEDTOOLS ]; then

    # get current version of bedtools
    VERSION=`bedtools --version | cut -f 2 -d '.'`

    # we want to have >= 27 in order to work correctly
    if [ "$VERSION" -lt "27"  ]; then
        install_bedtools
    fi
else
     install_bedtools
fi
