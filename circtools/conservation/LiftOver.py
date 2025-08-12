# This script is used to fetch the ortholog information per gene 
# using the REST API by ENSMBLE
import os, sys
import subprocess
import requests
import pybedtools
import platform

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
        
    def get_chain_file(self, from_id, to_id, dest_dir):
        tmp_name_species = to_id[0].upper() + to_id[1:]
        chain_filename = f"{from_id}To{tmp_name_species}.over.chain.gz"
        chain_path = os.path.join(dest_dir, chain_filename)

        if not os.path.exists(chain_path):
            # Try to download from UCSC
            url = f"https://hgdownload.cse.ucsc.edu/goldenPath/{from_id}/liftOver/{chain_filename}"
            print(f"Downloading chain file from {url} ...")
            response = requests.get(url, stream=True)
            if response.status_code == 200:
                os.makedirs(dest_dir, exist_ok=True)
                with open(chain_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"Downloaded chain file to {chain_path}")
            else:
                raise FileNotFoundError(f"Could not download chain file: {url}")

        return chain_path


    def call_liftover_binary(self):
    # Determine OS and architecture
        system = platform.system()
        machine = platform.machine()

        print(f"[DEBUG] Detected platform: {system} / {machine}", flush=True)

        # Identify correct liftOver subfolder
        if system == "Darwin":
            if machine == "x86_64":
                platform_dir = "AMD64/mac"
            elif machine == "arm64":
                platform_dir = "ARM64/mac"
            else:
                raise RuntimeError(f"Unsupported Mac architecture: {machine}")
        elif system == "Linux":
            platform_dir = "AMD64/linux"
        else:
            raise RuntimeError(f"Unsupported operating system: {system}")

        # Construct path to liftOver binary
        script_dir = os.path.dirname(os.path.realpath(__file__))
        parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
        liftover_utility = os.path.join(parent_dir, "contrib", "liftOver", platform_dir, "liftOver")
        print(f"[DEBUG] liftOver binary path: {liftover_utility}", flush=True)

        if not os.path.isfile(liftover_utility):
            raise FileNotFoundError(f"liftOver binary not found at: {liftover_utility}")

        # Build command
        command = [
            liftover_utility,
            self.liftover_input_file,
            self.chain_file,
            self.liftover_output_file,
            self.liftover_unlifted_file,
            "-multiple",
            "-minMatch=0.1"
        ]
        print(f"[DEBUG] Running liftOver command: {' '.join(command)}", flush=True)

        try:
            p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except Exception as e:
            print(f"[ERROR] Failed to start liftOver process: {e}", flush=True)
            raise

        try:
            out, err = p.communicate(timeout=60)  # 60s timeout to avoid infinite hang
        except subprocess.TimeoutExpired:
            p.kill()
            print(f"[ERROR] liftOver command timed out after 60s for {self.from_species} -> {self.to_species}", flush=True)
            return None

        p_status = p.wait()
        print(f"[DEBUG] liftOver exit code: {p_status}", flush=True)

        if out:
            print(f"[STDOUT]\n{out.decode(errors='replace')}", flush=True)
        if err:
            print(f"[STDERR]\n{err.decode(errors='replace')}", flush=True)

        if p_status != 0:
            print(f"[ERROR] liftOver failed with exit code {p_status}", flush=True)
            return None

        print(f"[DEBUG] liftOver completed successfully for {self.from_species} -> {self.to_species}", flush=True)
        return p


    
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
        self.chain_file = self.get_chain_file(self.from_id, self.to_id, os.path.join(self.tmpdir, "chain_files"))



        
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
        print(f"[DEBUG] Starting find_lifted_exons for {self.from_species} -> {self.to_species}", flush=True)

        lifted = self.parseLiftover()
        if lifted is None:
            print("[DEBUG] parseLiftover returned None", flush=True)
            return None

        chr = str(lifted[0])
        start = str(lifted[1])
        end = str(lifted[2])
        print(f"[DEBUG] Lifted coordinates: chr{chr}:{start}-{end}", flush=True)

        target_geneid = self.ortho_dict[self.to_species]
        print(f"[DEBUG] Target gene ID: {target_geneid}", flush=True)

        # First Ensembl request
        server = "https://rest.ensembl.org"
        ext = f"/overlap/region/{self.to_species}/{chr}:{start}-{end}?feature=gene;feature=exon"
        print(f"[DEBUG] Fetching exon overlaps from: {server+ext}", flush=True)

        try:
            r = requests.get(
                server+ext,
                headers={"Content-Type": "text/x-gff3"},
                timeout=(10, 60)  # 10s connect, 60s read
            )
        except requests.exceptions.RequestException as e:
            print(f"[ERROR] First exon fetch failed: {e}", flush=True)
            return None

        print(f"[DEBUG] First exon fetch HTTP status: {r.status_code}", flush=True)
        if not r.ok:
            print(f"[ERROR] Ensembl returned error {r.status_code}", flush=True)
            r.raise_for_status()
            sys.exit(1)

        print("WARNING! " + r.headers.get("X-RateLimit-Remaining", "?") + " REST API requests remaining!", flush=True)

        lifted_exons = self.parse_gff_rest(r.text)
        print(f"[DEBUG] Parsed {len(lifted_exons)} lifted exons", flush=True)

        lifted_exons_string = "\n".join(["\t".join(i) for i in lifted_exons])
        try:
            print("[DEBUG] Creating exon BedTool...", flush=True)
            exon_bed = pybedtools.BedTool(lifted_exons_string, from_string=True)
            region = "\t".join([chr, start, end])
            region_bed = pybedtools.BedTool(region, from_string=True)
            print("[DEBUG] Running exon_bed.intersect(region_bed)...", flush=True)
            intersect_exon = exon_bed.intersect(region_bed, wao=True)
            print("[DEBUG] Intersect completed", flush=True)
        except Exception as e:
            print(f"[ERROR] pybedtools intersect failed: {e}", flush=True)
            return None

        ortho_gene = self.ortho_dict[self.to_species]

        if str(intersect_exon).strip():
            intersect_out = [i.split("\t") for i in str(intersect_exon).strip().split("\n")]
            print(f"[DEBUG] Found {len(intersect_out)} intersecting exon entries", flush=True)
            intersect_out = [list(map(int, i)) for i in intersect_out]
            final_exon = sorted(intersect_out, key=lambda x: x[6], reverse=True)[0][:3]
            print(f"[DEBUG] Selected final exon: {final_exon}", flush=True)
            return final_exon
        else:
            print("[DEBUG] No intersecting exons found, fetching all exons for ortholog...", flush=True)
            ext = f"/overlap/id/{ortho_gene}?feature=exon"
            print(f"[DEBUG] Fetching all exons from: {server+ext}", flush=True)
            try:
                r = requests.get(
                    server+ext,
                    headers={"Content-Type": "text/x-gff3"},
                    timeout=(10, 60)
                )
            except requests.exceptions.RequestException as e:
                print(f"[ERROR] All exon fetch failed: {e}", flush=True)
                return None

            print(f"[DEBUG] All exon fetch HTTP status: {r.status_code}", flush=True)
            if not r.ok:
                r.raise_for_status()
                sys.exit(1)
            print("WARNING! " + r.headers.get("X-RateLimit-Remaining", "?") + " REST API requests remaining!", flush=True)

            all_exons = self.parse_gff_rest(r.text)
            print(f"[DEBUG] Parsed {len(all_exons)} total exons", flush=True)
            all_exons.sort(key=lambda x: x[1])

            try:
                all_exons_string = "\n".join(["\t".join(i) for i in all_exons])
                exon_bed_all = pybedtools.BedTool(all_exons_string, from_string=True)
                region = "\t".join([chr, start, end])
                region_bed = pybedtools.BedTool(region, from_string=True)
                print("[DEBUG] Running closest exon search...", flush=True)
                closest_exon = region_bed.closest(exon_bed_all, sortout=True)
                print("[DEBUG] Closest exon search completed", flush=True)
            except Exception as e:
                print(f"[ERROR] pybedtools closest failed: {e}", flush=True)
                return None

            final_exon = str(closest_exon).strip().split("\t")[-3:]
            print(f"[DEBUG] Selected closest exon: {final_exon}", flush=True)
            return final_exon


if __name__ == "__main__":
    lifted = liftover("human", "dog", ['2', '106145189', '106145475', 'UXS1', '0', '-'], "/scratch/circtools2/circtools/sample_data/temp", "test", 
                    {'dog': 'ENSCAFG00845009273', 'human': 'ENSG00000115652'}, "other")
    first_exon_liftover = lifted.find_lifted_exons()
