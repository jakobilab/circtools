Conservation module
********************************************************

The circtools padlock module is a specialized primer design tool tailored specifically for next-generation spatial transcriptomics technology. 

``circtools padlock`` is able to design padlock probes in batches of hundreds of circRNAs and linear RNAs based on circRNAs detected with ``circtools detect``, but can also work on lists of linear RNAs or circRNA isoforms or even entirely without any preliminary data purely based on the FASTA sequence of the circRNA.


General usage
--------------

A call to ``circtools conservation --help`` shows all available command line flags:

.. code-block:: bash

    usage: circtools [-h] -d DCC_FILE -g GTF_FILE -f FASTA_FILE [-O {mm,rn,hs,ss,cl}]
                     [-TS TARGET_SPECIES] [-s SEQUENCE_FILE] [-o OUTPUT_DIR] [-T EXPERIMENT_TITLE]
                     [-t GLOBAL_TEMP_DIR] [-G GENE_LIST [GENE_LIST ...]]
                     [-GL GENE_LIST_FILE [GENE_LIST_FILE ...]]
                     [-i ID_LIST [ID_LIST ...]] [-hg19] [-mm10] [-pairwise_flag]
    
    circular RNA conservation analysis
    
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
      -O {mm,rn,hs,ss,cl}, --organism {mm,rn,hs,ss,cl}
                            Organism of the study, mm =
                            Mus musculus, hs = Homo sapiens, rn = Rattus norvegicus,
                            ss = Sus scrofa, cl = Canis lupus familiaris
      -TS TARGET_SPECIES,	--target_species
                            List of target species IDs for which conservation score
                            needs to be calculated
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
      -i ID_LIST [ID_LIST ...], --id-list ID_LIST [ID_LIST ...]
                            Space-separated list of circRNA IDs. E.g. -i
                            "CAMSAP1_9_135850137_135850461_-"
                            "CAMSAP1_9_135881633_135883078_-"
      -hg19, --hg19
                            Are given circular co-ordinates for human from hg19 assembly?
                            If the flag is on, these will be converted into hg38.
      -mm10, --mm10
                            Are given circular co-ordinates for mouse from mm10 assembly?
                            If the flag is on, these will be converted into mm39.
      -pairwise_flag, --pairwise_flag
                            Should pairwise alignments be performed as well? 
                            Additional barplot will be plotted in this case.


Sample call to conservation module
----------------------------
A sample call to conservation using the `Jakobi et al. 2016 <https://www.sciencedirect.com/science/article/pii/S167202291630033X>`_ data requires the GTF file for exon information and the Fasta sequence of the reference genome in order to obtain the exon sequences.

.. code-block:: bash

    # run circtools conservation for circular RNA Slc8a1 to check its conservation in human and dog
    circtools conservation -d CircCoordinate -f Mus_musculus.GRCm38.dna.primary_assembly.fa -g Mus_musculus.GRCm38.90.gtf -O mm -G Slc8a1 -o test/ -t temp/ -TS hs -pairwise


.. code-block:: bash

	Start parsing GTF file
	Start merging GTF file outside the function
	Slc8a1_17_81647809_81649638_-
	extracting flanking exons for circRNA # 0 Slc8a1_17_81647809_81649638_-
	WARNING! 54986 REST API requests remaining!
	Processing target species:  hs
	*** Lifting over BSJ exon ***
	Successfully ran liftOver command human
	WARNING! 54985 REST API requests remaining!
	No nearby exon found. Trying for neaby exon search using orthology information.
	WARNING! 54984 REST API requests remaining!
	WARNING! 54983 REST API requests remaining!
	Lifted circle in target species  hs  is  ['2', '40097269', '40115629']
	mm(17:81647809- 0.000000
	hs(2:40097269-4 0.936299    0.000000
	    mm(17:81647809- hs(2:40097269-4
	Cleaning up




``circtools conservation`` takes a few seconds to process the input data. It fetches the information like gene orthologs, liftOver co-ordinates, exon sequences from REST API. The lifted over co-ordinates in target species are written in BED file. A phylogenetic tree for sequence alignement is drawn and saved in an SVG file.

If user wants to perform circle conservation analysis for species other than mentioned in the ``-O`` option, it can be easily done by editing the config file. An example config file is provided in the folder ``config/``. Following entries per species are required in order to include a new species:


.. code-block:: config

	mm:
  		input: # two letter abbreviation of the species (mm)
  		id: # genome versions for liftOver chain files (mm39)
 		name: # species alias according to Ensembl Rest API format (mouse)
  		ortho_id: # species name according to Ensembl Rest API format (mus_musculus)



