Installation
********************************************************


**circtools** is written in Python 3 (>=3.7) for most of the data processing code and R (>=4.0.0) for plotting and statistical analyses. The tool has a number of external dependencies, mostly standard bioinformatics tools and packages. The installation will, by default, try to install all required dependencies.

Installation is performed via `pip install circtools` or `python3 setup.py install`. No sudo access is required if the installation is suffixed with ``--user`` which will install the package in a user-writeable folder. In this case, the binaries should be installed to ``/home/$USER/.local/bin/`` (for Debian-based systems).


Supported operating systems
-----------------------------------

``circtools`` was developed and tested on Debian Buster (10) and Ubuntu Jammy Jellyfish (22.04). macOS is supported, but still in development. However, macOS functionality cannot be fully guaranteed yet.

Installation via PyPi
-----------------------------------

The default installation will install everything needed to run circtools *except R, STAR, or Stringtie* (see below). If you like you may install circtools locally (first call) or globally (second call, SU required).

.. code-block:: bash

    pip3 install circtools --user

Please note:

* The required R libraries will be installed in the default location in your home directory - unless you set enviromnet variable $R_LIBS_USER.
* In case want to install globally or into a dedicated 'venv' drop the --user option.

Installation via docker
-----------------------------------

Docker can be used as a simple alternative to install **circtools**, taking care of all Python and R dependencies, thus making it well-suited for novice users. The only requirement is a working docker installation. The following command will install the latest stable version of **circtools**:

.. code-block:: bash

    docker pull jakobilab/circtools

Subsequently, **circtools** can be run via

.. code-block:: bash

    docker run circtools [insert command args here]

.. code-block:: bash

Other than the additional command `docker run circtools` vs. just `circtools` all commands remain unchanged.


Installation via GitHub
--------------------------

The GitHub installation will install the most recent version directly from the source repository. Use this method if you want the latest fixes and features.

.. code-block:: bash

    git clone https://github.com/jakobilab/circtools.git
    cd circtools
    pip3 install . --verbose --user

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

Required dependencies
---------------------

External tools
^^^^^^^^^^^^^^^

* `bedtools [>= 2.27.1] <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_ required by the enrichment module

* `R [>= 4.0.0] <https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2>`_ required by visualisation scripts and the primer design module

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

Detailed manual installation
----------------------

Getting the source code
^^^^^^^^^^^^^^^^^^^^^^^

**Step 1**: Clone source code from GitHub:

.. code-block:: bash

    git clone https://github.com/jakobilab/circtools.git

Installation
^^^^^^^^^^^^

**Step 2**: Install circtools using the provided installation script. The ``--user`` flag installs circtools in your home folder, thus making sure you do not require any administrative rights during the installation:

.. code-block:: bash

    cd circtools
    pip3 install . install --verbose --user

R environment
^^^^^^^^^^^^^^

**Step 3**: Setting up R environment. In order for the automatic installation of R packages to work we need to set the package directory to a user-writeable path. The setup automatically sets that path to ``/home/$USER/.R/``.


Finishing up
^^^^^^^^^^^^

**Step 5**: Adding installation folder to ``$PATH``. In order for circtools to find all executables, the setup will add the folder ``/home/$USER/.local/bin/`` automatically to your ``.bashrc`` file

This closes the circtools installation. To verify that circtools has been correctly installed, try to call circtools for the first time:

.. code-block:: bash

    $> circtools --help
    usage: circtools [-V] <command> [<args>]
