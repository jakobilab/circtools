Padlock probe design module
********************************************************

The circtools padlock module is a specialized primer design tool tailored specifically for next-generation spatial transcriptomics technology. 

``circtools padlock`` is able to design padlock probes in batches of hundreds of circRNAs and linear RNAs based on circRNAs detected with ``circtools detect``, but can also work on lists of linear RNAs or circRNA isoforms or even entirely without any preliminary data purely based on the FASTA sequence of the circRNA.

Required tools and packages
----------------------------

``circtools primex`` depends on R, several R packages, and BioPython:

R packages:

* primex
* formattable
* kableExtra
* dplyr
* RColorBrewer
* colortools

Python libraries:

* BioPython>=1.71

All R package as well as Python dependencies are installed during the circtools installation.


General usage
--------------

A call to ``circtools padlock --help`` shows all available command line flags:

.. code-block:: bash

    usage: circtools [-h] -d DCC_FILE -g GTF_FILE -f FASTA_FILE [-O {mm,hs}]
                     [-s SEQUENCE_FILE] [-o OUTPUT_DIR] [-T EXPERIMENT_TITLE]
                     [-t GLOBAL_TEMP_DIR] [-G GENE_LIST [GENE_LIST ...]]
                     [-p PRODUCT_SIZE [PRODUCT_SIZE ...]]
                     [-i ID_LIST [ID_LIST ...]] [-j {r,n,f}] [-b]
    
    circular RNA primer design
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Input:
      -d DCC_FILE, --dcc-file DCC_FILE
                            CircCoordinates file from DCC / detect module
      -g GTF_FILE, --gtf-file GTF_FILE
                            GTF file of genome annotation e.g. ENSEMBL
      -f FASTA_FILE, --fasta FASTA_FILE
                            FASTA file with genome sequence (must match
                            annotation)
      -O {mm,hs}, --organism {mm,hs}
                            Organism of the study (used for primer BLASTing), mm =
                            Mus musculus, hs = Homo sapiens
      -s SEQUENCE_FILE, --sequence SEQUENCE_FILE
                            FASTA file containing the circRNA sequence (exons and
                            introns)
    
    Output options:
      -o OUTPUT_DIR, --output OUTPUT_DIR
                            Output directory (must exist)
      -T EXPERIMENT_TITLE, --title EXPERIMENT_TITLE
                            Title of the experiment for HTML output and file name
    
    Additional options:
      -t GLOBAL_TEMP_DIR, --temp GLOBAL_TEMP_DIR
                            Temporary directory (must exist)
      -G GENE_LIST [GENE_LIST ...], --genes GENE_LIST [GENE_LIST ...]
                            Space-separated list of host gene names. Primers for
                            CircRNAs of those genes will be designed.E.g. -G
                            "CAMSAP1" "RYR2"
      -p PRODUCT_SIZE [PRODUCT_SIZE ...], --product-size PRODUCT_SIZE [PRODUCT_SIZE ...]
                            Space-separated range for the desired PCR product.
                            E.g. -p 80 160 [default]
      -i ID_LIST [ID_LIST ...], --id-list ID_LIST [ID_LIST ...]
                            Space-separated list of circRNA IDs. E.g. -i
                            "CAMSAP1_9_135850137_135850461_-"
                            "CAMSAP1_9_135881633_135883078_-"
      -j {r,n,f}, --junction {r,n,f}
                            Should the forward [f] or reverse [r] primer be
                            located on the BSJ? [Default: n]
      -b, --no-blast        Should primers be BLASTED? Even if selected yes here,
                            not more than 50 primers willbe sent to BLAST in any
                            case.
      -r, --rna-type        Flag for RNA type for which you want to                                                                                                                 generate padlock probes . 0 for Circular RNAs only, 1 for
                            Linear RNAs only and 2 for Both. DEFAULT 2
      -svg, --no-svg        Should the SVG files for graphical                                                         representation be generated?



