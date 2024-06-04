# This script is used to fetch the ortholog information per gene 
# using the REST API by ENSMBLE
import os, sys
import subprocess

class liftover(object):

    def __init__(self, from_species, to_species, bed_coord, tmpdir, prefix, flag) -> None:
        self.from_species = from_species
        self.to_species = to_species
        self.from_coord = bed_coord     # BED coordinates in form of a list of chr, start and stop, score and strand
        self.tmpdir = tmpdir
        self.prefix = prefix
        self.flag = flag

    def call_liftover_binary(self):
        # encapsulated liftover binary call

        # liftover command
        liftover_utility = "/home/skulkarni/liftOver"
        command = liftover_utility + " " + self.liftover_input_file + " " + self.chain_file + " " + self.liftover_output_file + " " + self.liftover_unlifted_file + "  -multiple -minMatch=0.1"
        print(command)
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        return(p)

    
    def lifting(self):
        # function to perform actual lifting

        ## check if the flag for mm10 and hg19 conversion is true. 
        if self.flag == "mm10":
            # this is only internal leftover for mouse from version mm10 to mm39
            self.from_id = "mm10"
            self.to_id = "mm39"
            tmp_from_bed = self.tmpdir + self.prefix + "_liftover_internal.tmp"
            open(tmp_from_bed, 'w').close()             # erase old contents
        
        elif self.flag == "hg19":
            # this is only internal leftover for mouse from version hg19 to hg38
            self.from_id = "hg19"
            self.to_id = "hg38"
            tmp_from_bed = self.tmpdir + self.prefix + "_liftover_internal.tmp"
            open(tmp_from_bed, 'w').close()             # erase old contents
        
        elif self.flag == "other":
            species_IDs_dict = {"mouse":"mm39", "human":"hg38", "pig":"susScr11", "dog":"canFam6"}
            self.from_id = species_IDs_dict[self.from_species]
            self.to_id = species_IDs_dict[self.to_species]
            tmp_from_bed = self.tmpdir + self.prefix + "_liftover.tmp"
            open(tmp_from_bed, 'w').close()             # erase old contents
        
        else:
            print("Unidentified flag for liftOver function:", self.flag)
            sys.exit()

        
        with open(tmp_from_bed, 'a') as data_store:
            data_store.write("chr" + "\t".join(self.from_coord) + "\n")
        # chain file
        chain_file = self.from_id + "To" + self.to_id.title() + ".over.chain.gz"
        self.chain_file = chain_file
        
        tmp_to_bed = tmp_from_bed + ".out"              # output file
        open(tmp_to_bed, 'a').close()                   # erase old contents
        tmp_unlifted = tmp_from_bed + ".unlifted"       # unlifted file
        open(tmp_unlifted, 'a').close()                   # erase old contents

        self.liftover_input_file = tmp_from_bed
        self.liftover_output_file = tmp_to_bed
        self.liftover_unlifted_file = tmp_unlifted

        # liftover binary call
        p = self.call_liftover_binary()
        
        # check the command status 
        out, err = p.communicate()
        p_status = p.wait()
        if (p_status != 0):
            print("liftOver command not run successfully. Exiting!")
            print(out, err)
            sys.exit()
        else:
            print("Successfully ran liftOver command " + self.to_species)

    def parseLiftover(self):
        # function to parse liftover output files and return the lifted coordinates to main function
        self.lifting()
        tmp_from_bed = self.liftover_input_file
        tmp_to_bed = self.liftover_output_file
        tmp_unlifted = self.liftover_unlifted_file
        
        # read in the unlifted file to see if there were any errors 
        # if not, read the output file and print the lifted coordinates
        if os.stat(tmp_unlifted).st_size != 0:
            print("Unlifted coordinates present. Liftover did not run well. Exiting!")
            sys.exit()
        else:
            fin = open(tmp_to_bed).readlines() #.strip().split("\t")
            with open(tmp_to_bed) as fin:
                lines = fin.read().splitlines()
            if (len(lines) == 1):
                #print(lines)
                lifted_coordinates = lines[0].split("\t")
            else:
                # somehow the lifted coordinates are split into two. 
                for line in lines:
                    print(line)
            
            lifted_coordinates[0] = lifted_coordinates[0].replace("chr", "")
            print("Lifted coordinates:", lifted_coordinates)
        return(lifted_coordinates)