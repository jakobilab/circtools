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
        import time
        print("[DEBUG] Starting fetch_info", flush=True)
        print(f"[DEBUG] from_species: {self.from_species}", flush=True)
        print(f"[DEBUG] gene_symbol: {self.gene_symbol}", flush=True)
        print(f"[DEBUG] dict_species_ortholog: {self.dict_species_ortholog}", flush=True)

        # define the species for which you are finding orthologs
        species = list(self.dict_species_ortholog.values())
        print(f"[DEBUG] species list: {species}", flush=True)

        to_species = list(set(species) - set([self.from_species]))
        print(f"[DEBUG] to_species (raw): {to_species}", flush=True)

        to_species = ["target_species=" + str(i) for i in to_species]
        str_to_species = ";".join(to_species)
        print(f"[DEBUG] str_to_species param: {str_to_species}", flush=True)

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
        print(f"[DEBUG] Final URL: {url}", flush=True)

        # Retry loop without importing Retry/HTTPAdapter
        max_attempts = 5
        backoff = 2
        r = None
        for attempt in range(1, max_attempts + 1):
            print(f"[DEBUG] Attempt {attempt} starting...", flush=True)
            start_time = time.time()
            try:
                r = requests.get(
                    url,
                    headers={"Content-Type": "application/json"},
                    timeout=(10, 60),  # 10s connect, 60s read
                    verify=True
                )
                elapsed = time.time() - start_time
                print(f"[DEBUG] Attempt {attempt} completed in {elapsed:.2f}s", flush=True)
                if r.ok:
                    print(f"[DEBUG] HTTP {r.status_code} OK", flush=True)
                    break
                else:
                    print(f"[DEBUG] HTTP {r.status_code} {r.reason}, retrying in {backoff}s...", flush=True)
            except requests.exceptions.RequestException as e:
                elapsed = time.time() - start_time
                print(f"[DEBUG] Attempt {attempt} failed after {elapsed:.2f}s: {type(e).__name__} -> {e}", flush=True)

            if attempt < max_attempts:
                time.sleep(backoff)
                backoff *= 2  # exponential backoff
        else:
            print("[DEBUG] Could not fetch ortholog information from ENSEMBL after retries. Exiting!", flush=True)
            sys.exit(1)

        remaining = r.headers.get("X-RateLimit-Remaining", "?")
        print(f"WARNING! {remaining} REST API requests remaining!", flush=True)
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
