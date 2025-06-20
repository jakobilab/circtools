# This script is used to fetch the ortholog information per gene 
# using the REST API by ENSMBLE

import sys
import requests

class fetch(object):

    def __init__(self, gene_symbol, from_species, dict_species_ortholog, target_species_code) -> None:
        self.gene_symbol = gene_symbol
        self.from_species = from_species
        self.dict_species_ortholog = dict_species_ortholog
        self.target_species_code = target_species_code  # e.g., 'hs'

    def fetch_info(self):
        # Map the short code (like 'hs') to the full species name (like 'human')
        if self.target_species_code not in self.dict_species_ortholog:
            raise ValueError(f"Target species code '{self.target_species_code}' not found in species mapping.")

        target_species = self.dict_species_ortholog[self.target_species_code]
        str_to_species = f"target_species={target_species}"

        server = "https://rest.ensembl.org"
        ext = f"/homology/symbol/{self.from_species}/{self.gene_symbol}?format=condensed;type=orthologues;{str_to_species}"

        try:
            r = requests.get(server + ext, headers={ "Content-Type": "application/json" })
        except requests.exceptions.RequestException as e:
            raise SystemExit(e)

        if not r.ok:
            r.raise_for_status()
            print("Could not fetch ortholog information from ENSEMBL. Exiting!")
            sys.exit()

        print("WARNING! " + r.headers.get("X-RateLimit-Remaining", "Unknown") + " REST API requests remaining!")
        return r

    def parse_json(self):
        ortho_dict = {}
        species_dict = self.dict_species_ortholog

        out_string = self.fetch_info().json()
        original_species_gene_id = out_string['data'][0]["id"]
        ortho_dict[self.from_species] = original_species_gene_id

        ortho_dir = out_string['data'][0]["homologies"]
        for each_key in ortho_dir:
            species_code = [k for k, v in species_dict.items() if v == each_key["species"]]
            if species_code:
                ortho_dict[species_code[0]] = each_key["id"]

        return ortho_dict
    
#obj = fetch("SLC8A1", "human")
#obj.fetch_info()
#obj.parse_json()
