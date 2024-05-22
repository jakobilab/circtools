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

        tmp_from_bed = self.tmpdir + self.prefix + "_circtools_flanking_exons.tmp"
        open(tmp_from_bed, 'w').close()             # erase old contents
        with open(tmp_from_bed, 'a') as data_store:
            data_store.write("\t".join(self.from_coord))
        # chain file
        chain_file = self.from_id + "To" + self.to_id + ".over.chain.gz"
        # output file
        tmp_to_bed = tmp_from_bed + ".out.tmp"
        # unlifted file
        tmp_unlifted = tmp_from_bed + ".unlifted.tmp" 

        liftover_utility = "/home/skulkarni/liftOver"
        command = liftover_utility + " " + chain_file + " " + tmp_to_bed + " " + tmp_unlifted 
        print(command)