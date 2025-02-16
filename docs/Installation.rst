Installation
********************************************************


**circtools** is written in Python 3 (>=3.8) for most of the data processing code and R (>=4.0.0) for plotting and statistical analyses. The tool has a number of external dependencies, mostly standard bioinformatics tools and packages. The installation will, by default, try to install all required dependencies.

Installation is performed via `pip install circtools` or `python3 setup.py install`. No sudo access is required if the installation is suffixed with ``--user`` which will install the package in a user-writeable folder. In this case, the binaries should be installed to ``/home/$USER/.local/bin/`` (for Debian-based systems).


Supported operating systems
-----------------------------------

``circtools`` was developed and tested on Debian Bookworm (12), Ubuntu Jammy Jellyfish (22.04), and Ubuntu Noble Numba (24.04). macOS is supported via docker installation. However, macOS functionality running the pip installation method cannot be fully guaranteed yet.

Installation via docker
-----------------------------------

The latest circtools docker version will be downloaded directly from GitHub. The container contains `all` dependencies required to run `circtools` except STAR and Bowtie.

.. code-block:: console

    docker pull ghcr.io/jakobilab/circtools/circtools:latest

An bash alias to call circtools "natively" and skip the unwieldy full docker command is recommended:

.. code-block:: console

    alias circtools='docker run --rm -v "`pwd`":/circtools/ ghcr.io/jakobilab/circtools/circtools'

This line can be added to the `.bashrc` or `.profile` file to be automatically loaded after login.


Installation via PyPi
-----------------------------------

The default installation will install everything needed to run circtools *except R, STAR, or Stringtie* (see below).  We recommend to install circtools in a virtual environment (venv) to separate dependencies from other packages as well as the OS.

.. code-block:: bash

    python3 -m venv circtools # create virtual environment
    source circtools/bin/activate # activate virtual environment
    python3 -m pip install circtools # install latest circtools version from pypy

.. note::

    The required R libraries will be installed in the default location in the home directory - unless the environment variable $R_LIBS_USER is set.


Installation via GitHub
--------------------------

The GitHub installation will install the most recent version directly from the source repository. Use this method if you want the latest fixes and features.

.. code-block:: bash

    git clone https://github.com/jakobilab/circtools.git
    cd circtools
    python3 -m venv circtools_venv # create virtual environment
    source circtools_venv/bin/activate # activate virtual environment
    python3 -m pip install . # install latest circtools version


Installation of R dependencies
--------------------------------

All R packages required by circtools can be automatically installed with a single command:

.. code-block:: bash

    circtools_install_R_dependencies

The R packages require certain development libraries installed in order to be compiled from source.

The following libraries are required (Ubuntu/Debian package names given):

- libpng-dev
- zlib1g-dev
- libbz2-dev
- libjpeg-turbo8-dev
- libcurl4-openssl-dev
- libxml2-dev
- libblas-dev
- liblzma-dev
- libfontconfig1-dev
- liblapack-dev
- libssl-dev
- libharfbuzz-dev
- libfribidi-dev
- libfreetype6-dev
- libtiff5-dev
- libjpeg-dev

A simple command to install all of these libraries on an Ubuntu/Debian system would be:

.. code-block:: bash

    apt-get install --no-install-recommends r-base python3 python3-dev make g++ gfortran libpng-dev zlib1g-dev libbz2-dev libjpeg-turbo8-dev libcurl4-openssl-dev libxml2-dev libblas-dev liblzma-dev libfontconfig1-dev liblapack-dev libssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libtiff5-dev libjpeg-dev

The command above only installs the minimal required packages, no other recommend packages are installed to keep the system lean.



Updating circtools
--------------------------

You may want to update the circtools package if new versions are published. Similar to the initial installation, there are two ways to update circtools:

.. code-block:: bash

    pip3 install circtools --user --upgrade

.. code-block:: bash

    cd /path/to/circtools/repo/
    git pull
    cd circtools/
    pip3 install . install --verbose --user --upgrade


Finishing up
------------------
In order for circtools to find all executables, the setup will add the folder ``/home/$USER/.local/bin/`` automatically to your ``.bashrc`` file

This closes the circtools installation. To verify that circtools has been correctly installed, try to call circtools for the first time:

.. code-block:: bash

    $> circtools --help
    usage: circtools [-V] <command> [<args>]


Required dependencies
---------------------

External tools
^^^^^^^^^^^^^^^

* `bedtools [>= 2.27.1] <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_ required by the enrichment module

* `R [>= 4.0.0] <https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-20.04>`_ required by visualisation scripts and the primer design module

* `STAR [>= 2.6.0] <https://github.com/alexdobin/STAR>`_ required by the ``detect`` and ``reconstruct`` module to map RNA-seq reads against a reference genome and detect back splice junctions

* `Stringtie [>= 1.3.3b, optional] <https://github.com/gpertea/stringtie>`_ required by the ``exon`` module to carry out exon level analyses.

The primer design module as well as the exon analysis and circRNA testing module require a working installation of `R <https://cran.r-project.org/>`_ with `BioConductor <https://www.bioconductor.org/install/>`_. All R packages required are automatically installed during the setup.

.. important:: The setup scripts assumes that the folder for R plugins is writeable (either in the user's home or the system folder).

Required Python packages (automatically installed)
^^^^^^^^^^^^^^^
- HTSeq >= 0.11.0
- pysam >= 0.16.0.1
- numpy >= 1.14.5
- pybedtools >= 0.7.10
- biopython >= 1.71
- scipy >= 0.19.0
- reportlab >= 3.3.0
- pandas >= 0.25.0
- statsmodels >= 0.9.0
