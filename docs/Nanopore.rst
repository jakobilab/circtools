Nanopore module
********************************************************

The circtools nanopore module is designed to process



``circtools nanopore`` requires sequencing reads that have been produced using the protocol outlined in `Rahimi et al. (2021) <https://doi.org/10.1038/s41467-021-24975-z>`_. Specifically, OxFord Nanopore sequencing reds from this protocol will cover circular RNAs over their full length once or even multiple times. The circtools ``circtools nanopore`` module is based on code published as part of the `Rahimi et al. (2021) <https://doi.org/10.1038/s41467-021-24975-z>`_ study (`GitHub repository <https://github.com/omiics-dk/long_read_circRNA>`_). The original code has been thoroughly updated and the ``circtools nanopore`` module now provides support for `hg38` and `mm10` through a completely new reference data retrieval system that automatically downloads and postprocesses reference data from public databases such as Genome Browser and ENSEMBL without the need for short-term links to private online drives. This new system also enables easy addition of other model species as long as they are supported by ENSEMBL and the Genome Browser. Moreover, Perl scripts were replaced by new implementations in Python3 to seamlessly integrate in the circtools 2.0 software framework and reduce the overall software maintenance burden.


Required tools and packages
--------------------------------
``nanopore`` depends several external tools, namely

* `Nanofilt <https://github.com/wdecoster/nanofilt/releases>`_: Quality filtering of Nanopore data
* `bedtools <https://github.com/arq5x/bedtools2/releases>`_: Genome arithmetics
* `pblat <https://github.com/icebert/pblat/releases>`_: Alignment of Nanopore reads
* `samtools <https://github.com/samtools/samtools/releases>`_: Alignment manipulation


General usage
--------------

A call to ``circtools nanopore --help`` shows all available command line flags:

.. code-block:: bash

    usage: circtools.py [-h] (-r | -c | -d) [-s SAMPLE] [-R REFERENCE_PATH] [-O OUTPUT_PATH] [-C {hg19,hg38,mm9,mm10}] [-t THREADS] [-D] [-k]

    circular RNA detection in Oxford Nanopore data

    options:
      -h, --help            show this help message and exit
      -r, --run             Run the analysis
      -c, --check           Check the installation for required software.
      -d, --download        Download third-party data, such as genomes required for the analysis.

    Options:
      -s SAMPLE, --sample SAMPLE
                            Provide a sample input .fq.gz file that should be processed.
      -R REFERENCE_PATH, --reference-path REFERENCE_PATH
                            Provide a path for where the reference data is located. Default is './data'.
      -O OUTPUT_PATH, --output OUTPUT_PATH
                            Provide a path for where the output data is stored.
      -C {hg19,hg38,mm9,mm10}, --config {hg19,hg38,mm9,mm10}
                            Required. Select which genome build the sample that is from, and specify which genome reference files should be used.
      -t THREADS, --threads THREADS
                            Number of threads for parallel steps. Default: 4.
      -D, --dry-run         Perform all of the input checks without starting the detection scripts.
      -k, --keep-temp       Keep all of the temporary files.




Setup: Check if external software is available
^^^^^^^^^^^^

.. note::

    If the Docker image is used all required software is already installed within the image.


In order to check if the external software has been installed correctly and can be used by the nanopore module a check can be run:

.. code-block:: bash

    circtools nanopore -c

This should produce the following output, indicating tat software is accessible:

.. code-block:: bash

    Checking for bedtools
    Checking for NanoFilt
    Checking for pblat
    Checking for samtools

    All of the expected software requirements are present!

Should software not be installed, e.g. pblat, an error message is shown:

.. code-block:: bash

    Checking for bedtools
    Checking for NanoFilt
    Checking for pblat
            Unable to find pblat!
    Checking for samtools

    ERROR: Some of the required software is missing!


Step 1: Download required data
^^^^^^^^^^^^
.. code-block:: bash

    circtools nanopore -d -R reference/ -C hg38

Here the reference data will be downloaded into in the folder ``reference/`` and we are download all require files for the human genome, build `hg38`. The folder will be automatically created if it does not exist. For each reference genome, a suitable sub-folder will be created, e.g. `hg38` which contains all required and post-processed files. All downloads are linking to public sources, such as the `Genome Browser <https://genome.ucsc.edu/>`_; links are stored in YAML files available in the `GitHub repository <https://github.com/jakobilab/circtools/tree/master/circtools/nanopore/config>`_.

We are welcoming pull requests for additional genome builds!

The download progress is visible in the command line together with automatic post-processing:

.. code-block:: bash

    Storing reference data in reference/
    Downloading genome.fa.gz: 100%|█████████████████████████████████████████████| 984M/984M [01:00<00:00, 16.3MB/s]
    Unpacking.
    Done.
    Downloading genome.chrom.sizes: 100%|██████████████████████████████████████| 11.7k/11.7k [00:00<00:00, 602kB/s]
    Downloading refFlat.csv.gz: 3.92MB [00:01, 3.22MB/s]
    Creating refFlat-based exon files
    Downloading gencode.csv.gz: 100%|█████████████████████████████████████████| 59.0M/59.0M [00:34<00:00, 1.70MB/s]
    Unpacking.
    Done.
    Creating GENCODE-based exon files
    Start parsing GTF file
    Downloading gencode_intron.bed.gz: 8.74MB [00:03, 2.51MB/s]
    Unpacking.
    Done.
    Downloading est.bed.gz: 444MB [03:20, 2.22MB/s]
    Unpacking.
    Done.

In the above example, the folder ``reference/hg38/`` should now contain the following files occupying around 8GB of disk space.

.. code-block:: bash

    est.bed
    gencode.csv
    gencode.csv.exon.bed
    gencode.csv.exon.merge.bed
    gencode_intron.bed
    genome.chrom.sizes
    genome.fa
    refFlat.csv.gz
    refFlat.csv.merged.bed
    refFlat.csv.sort.bed
    refFlat.csv.unique.bed

The file names are identical for each genome build, only the folder name indicates which genome is stored in each folder.

Step 2: Run the nanopore pipeline
^^^^^^^^^^^^

To run the main workflow of the ``circtools nanopore`` module, users need to specify the reference genome (``-R reference/``), output path (``-O results/``), and the FASTQ file containing the Oxford Nanopore reads (``-s human_nanopore.fastq.gz``). An example dataset consisting of `100k human brain nanopore reads is available for download <https://github.com/jakobilab/circtools/raw/refs/heads/master/tests/data/human_nanopore.fastq.gz>`_. The ``--threads 16`` argument is optional, but can be supplied to speed up processing by using multiple CPU threads, in this case 16 threads:

.. code-block:: bash

    circtools.py nanopore -r -s human_nanopore.fastq.gz -R reference/ -C hg38 -O results/ --threads 16

The pipeline outputs a number of output files, specifically:

.. code-block:: bash

    ls -la results/

    human_nanopore.circ_circRNA_exon_usage_length_of_exons.txt
    human_nanopore.circRNA_candidates.annotated.txt
    human_nanopore.novel.cryptic.spliced.exons.txt
    human_nanopore.novel.exons.2reads.filter.bed
    human_nanopore.novel.exons.2reads.phases.tab
    human_nanopore.Potential_multi-round_circRNA.fa
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.10reads.bed
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.20reads.bed
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.2reads.bed
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.3reads.bed
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.50reads.bed
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.5reads.bed
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.bed
    human_nanopore.scan.circRNA.psl.split.merge.flank2.allExons.notGencode.bed
    human_nanopore.scan.Potential_multi-round_circRNA.psl.annot.bed
    human_nanopore.scan.Potential_multi-round_circRNA.psl.annot.count.txt

The files are prefixed with the sample name (input FASTQ file name minus extension) and are named intuitively. The main output file has the suffix `circRNA_candidates.annotated.txt` and contains the list of circRNAs detected in the run. Specifically, the files contains the following columns for each circRNA:

.. code-block:: bash

     1  internal_circRNA_name
     2  chr
     3  start
     4  end
     5  description
     6  BSJ_reads
     7  strand
     8  gene
     9  reserved
    10  reserved
    11  reserved
    12  mean_read_coverage
    13  mean_gene_coverage
    14  mean_exon_coverage
    15  mean_EST_coverage
    16  mean_intron_coverage
    17  min_exon_adjust
    18  max_exon_adjust
    19  mean_exon_adjust

