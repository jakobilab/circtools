# This script is used to fetch the ortholog information per gene 
# using the REST API by ENSMBLE
import os, sys
import subprocess
import requests
import pybedtools

class liftover(object):

    def __init__(self, from_species, to_species, bed_coord, tmpdir, prefix, orthologs, flag, dict_species_liftover) -> None:
        self.from_species = from_species
        self.to_species = to_species
        self.from_coord = bed_coord     # BED coordinates in form of a list of chr, start and stop, score and strand
        self.gene_name = bed_coord[3]
        self.tmpdir = tmpdir
        self.prefix = prefix
        self.flag = flag
        self.ortho_dict = orthologs
        self.dict_species_liftover = dict_species_liftover

    def call_liftover_binary(self):
        # encapsulated liftover binary call

        # liftover command
        liftover_utility = os.path.dirname(os.path.realpath(__file__)) + "/liftOver"
        command = liftover_utility + " " + self.liftover_input_file + " " + self.chain_file + " " + self.liftover_output_file + " " + self.liftover_unlifted_file + "  -multiple -minMatch=0.1"
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
            #species_IDs_dict = {"mouse":"mm39", "human":"hg38", "pig":"susScr11", "dog":"canFam6", "rat":"rn7"}
            species_IDs_dict = self.dict_species_liftover
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
        tmp_name_species = self.to_id[0].upper() + self.to_id[1:]
        chain_file = self.from_id + "To" + tmp_name_species + ".over.chain.gz"
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
            #sys.exit()
            return(None)
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
                    print("Lifted coordinates are splitted into two regions", line)
            
            lifted_coordinates[0] = lifted_coordinates[0].replace("chr", "")
            #print("Lifted coordinates:", lifted_coordinates)
            return(lifted_coordinates)
    
    def parse_gff_rest(self, output):
        # function to parse the gff output from REST API exon extraction information
        # returns list of fetched exons overlapping a given region of interest
        exon_list = []
        out = output.strip().split("\n")
        for eachline in out:
            if (eachline.startswith("#")):  continue
            eachline = eachline.split("\t")
            if (eachline[2] == "exon"):
                # reduce 1bp from start because circcoordinates are 0 based
                eachline[3] = str(int(eachline[3]) -1 )
                exon_list.append([eachline[0], eachline[3], eachline[4]])
        
        exon_list = [list(x) for x in set(tuple(x) for x in exon_list)] 
        return(exon_list)

    def find_lifted_exons(self):
        # function to see if the lifted coordinates are exons. If not, take nearby exons

        lifted = self.parseLiftover()
        if lifted == None:
            return(None)
        chr = str(lifted[0])
        start = str(lifted[1])
        end = str(lifted[2])

        #print("To be lifted coordinates: ", self.from_coord)

        target_geneid = self.ortho_dict[self.to_species]
        #print("Extracting exons for : ", target_geneid)
        server = "https://rest.ensembl.org"
        ext = "/overlap/region/" + self.to_species + "/" + chr + ":" + start + "-" + end + "?feature=gene;feature=exon"

        try:
            r = requests.get(server+ext, headers={ "Content-Type" : "text/x-gff3"})
        except requests.exceptions.RequestException as e:
            raise SystemExit(e)
        
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        print("WARNING! "+ r.headers["X-RateLimit-Remaining"] + " REST API requests remaining!")
        
        # call above gff parsing function on this output
        lifted_exons = self.parse_gff_rest(r.text)
        #print("Lifted exons:", lifted_exons)

        # now perform bedtools operation to find out the correct exon boundaries
        lifted_exons_string = "\n".join(["\t".join(i) for i in lifted_exons])
        exon_bed = pybedtools.BedTool(lifted_exons_string, from_string = True)
        region = "\t".join([chr, start, end])
        region_bed = pybedtools.BedTool(region, from_string = True)
        intersect_exon = exon_bed.intersect(region_bed, wao=True)
        #print("Intersect:", str(intersect_exon))
        ortho_gene = self.ortho_dict[self.to_species]

        if (intersect_exon != ""):
            # parse the above information to keep the longest overlapping exon
            intersect_out = [i.split("\t") for i in str(intersect_exon).strip().split("\n")]
            intersect_out = [list(map(int, i)) for i in intersect_out]
            #print(intersect_out)

            # sort the above list based on 7th element i.e. overlap bases and take the exon corresponding to the maximum overlap
            final_exon = sorted(intersect_out, key=lambda x: x[6], reverse=True)[0][:3]
            #print("Final:", final_exon)

            return(final_exon)                  # the sequences will be extracted for this exon
        else:
            # no intersecting exons found. In this case, fetch orthologs for this species and extract exon information
            # and intersect with exons to take the closest exon information
            print("No nearby exon found. Trying for neaby exon search using orthology information.")

            # fetch the exon information for this genee

            server = "https://rest.ensembl.org"
            ext = "/overlap/id/" + ortho_gene + "?feature=exon"
            #print(server+ext)
            try:
                r = requests.get(server+ext, headers={ "Content-Type" : "text/x-gff3"})
            except requests.exceptions.RequestException as e:
                raise SystemExit(e)
            
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            print("WARNING! "+ r.headers["X-RateLimit-Remaining"] + " REST API requests remaining!")
            #print("All exons for gene id "+ ortho_gene)
            #print(r.text)
            
            # call above gff parsing function on this output
            all_exons = self.parse_gff_rest(r.text)
            all_exons.sort(key=lambda x: x[1])

            # now perform bedtools operation to find out the correct exon boundaries
            all_exons_string = "\n".join(["\t".join(i) for i in all_exons])
            #print("All exons for closest:", all_exons_string)
            exon_bed_all = pybedtools.BedTool(all_exons_string, from_string = True)
            region = "\t".join([chr, start, end])
            #print("To  find closest region from :", region)
            region_bed = pybedtools.BedTool(region, from_string = True)
            closest_exon = region_bed.closest(exon_bed_all, sortout=True)
            final_exon = str(closest_exon).strip().split("\t")[-3:]
            #print("Closest exon:",   final_exon)

            return(final_exon)

if __name__ == "__main__":
    lifted = liftover("human", "dog", ['2', '106145189', '106145475', 'UXS1', '0', '-'], "/scratch/circtools2/circtools/sample_data/temp", "test", 
                    {'dog': 'ENSCAFG00845009273', 'human': 'ENSG00000115652'}, "other")
    first_exon_liftover = lifted.find_lifted_exons()
