# This script is used to fetch the ortholog information per gene 
# using the REST API by ENSMBLE

class liftover(object):

    def __init__(self, from_species, to_species, bed_coord, tmpdir, prefix) -> None:
        self.from_species = from_species
        self.to_species = to_species
        self.from_coord = bed_coord     # BED coordinates in form of a list of chr, start and stop, score and strand
        self.tmpdir = tmpdir
        self.prefix = prefix

    def lifting(self):
        # function to perform actual lifting

        species_IDs_dict = {"mouse":"mm10", "human":"hg38", "pig":"susScr11", "dog":"canFam6"}
        self.from_id = species_IDs_dict[self.from_species]
        self.to_id = species_IDs_dict[self.to_species]

        tmp_from_bed = self.tmpdir + self.prefix + "_liftover.tmp"
        open(tmp_from_bed, 'w').close()             # erase old contents
        with open(tmp_from_bed, 'a') as data_store:
            data_store.write("chr" + "\t".join(self.from_coord) + "\n")
        # chain file
        chain_file = self.from_id + "To" + self.to_id.title() + ".over.chain.gz"
        
        tmp_to_bed = tmp_from_bed + ".out"              # output file
        open(tmp_to_bed, 'a').close()                   # erase old contents
        tmp_unlifted = tmp_from_bed + ".unlifted"       # unlifted file
        open(tmp_unlifted, 'a').close()                   # erase old contents

        liftover_utility = "/home/skulkarni/liftOver"
        command = liftover_utility + " " + tmp_from_bed + " " + chain_file + " " + tmp_to_bed + " " + tmp_unlifted + "  -multiple -minMatch=0.1"
        print(command)
