# This script is used to fetch the ortholog information per gene 
# using the REST API by ENSMBLE

import sys
import requests

class fetch(object):

    def __init__(self, gene_symbol, from_species, dict_species_ortholog) -> None:
        self.gene_symbol = gene_symbol
        self.from_species = from_species
        self.dict_species_ortholog = dict_species_ortholog

    def fetch_info(self):
    # define the species for which you are finding orthologs
        species = list(self.dict_species_ortholog.values())
        to_species = list(set(species) - set([self.from_species]))
        to_species = ["target_species=" + str(i) for i in to_species]
        str_to_species = ";".join(to_species)

        server = "https://rest.ensembl.org"
        ext = (
            "/homology/symbol/"
            + self.from_species
            + "/"
            + self.gene_symbol
            + "?format=condensed;type=orthologues;"
            + str_to_species
        )
        url = server + ext

        # Retry loop without importing Retry/HTTPAdapter
        max_attempts = 5
        backoff = 2
        for attempt in range(1, max_attempts + 1):
            try:
                r = requests.get(
                    url,
                    headers={"Content-Type": "application/json"},
                    timeout=(10, 60)  # 10s connect, 60s read
                )
                if r.ok:
                    break
                else:
                    print(f"Attempt {attempt}: HTTP {r.status_code}, retrying in {backoff}s...")
            except requests.exceptions.RequestException as e:
                print(f"Attempt {attempt} failed: {e}, retrying in {backoff}s...")
            if attempt < max_attempts:
                import time
                time.sleep(backoff)
                backoff *= 2  # exponential backoff
        else:
            print("Could not fetch ortholog information from ENSEMBL after retries. Exiting!")
            sys.exit(1)

        remaining = r.headers.get("X-RateLimit-Remaining", "?")
        print(f"WARNING! {remaining} REST API requests remaining!")
        return r


    def parse_json(self):
        ortho_dict = {}
        # species_dict = {"canis_lupus_familiaris": "dog", "mus_musculus": "mouse", "homo_sapiens": "human",
        #                "sus_scrofa": "pig", "rattus_norvegicus": "rat"}
        species_dict = self.dict_species_ortholog
        # this function takes JSON output from REST API and parses the ortholog information
        out_string = self.fetch_info().json()
        original_species_gene_id = out_string['data'][0]["id"]
        ortho_dict[self.from_species] = original_species_gene_id 
        ortho_dir = out_string['data'][0]["homologies"]
        for each_key in ortho_dir:
            ortho_dict[species_dict[each_key["species"]]] = each_key["id"] 
        
        return ortho_dict
    
#obj = fetch("SLC8A1", "human")
#obj.fetch_info()
#obj.parse_json()
