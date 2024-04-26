# This script is used to fetch the ortholog information per gene 
# using the REST API by ENSMBLE

import os, sys
import requests

class fetch(object):
    def __init__(self) -> None:
        pass

    def fetch_info(self, gene_symbol, from_species):
        # define the species for which you are finding orthologs
        species = ["mouse", "human", "rat", "pig", "dog"]
        to_species = list(set(species)-set([from_species]))
        to_species = ["target_species="+str(i) for i in to_species]
        str_to_species = ";".join(to_species)

        server = "https://rest.ensembl.org"
        ext = "/homology/symbol/" + from_species + "/" + gene_symbol + "?format=condensed;type=orthologues;" + str_to_species
 
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
        if not r.ok:
            r.raise_for_status()
            sys.exit()
  
        print(r.text)

    def parse_json(self, r):
        # this function takes JSON output from REST API and fetches the orthologs information
        out_string = r.json()
        original_species_id = out_string['data'][0]["id"]
        ortho_dir = out_string['data'][0]["homologies"]
        for each_key in ortho_dir.keys():
            print(each_key, )

obj = fetch()
obj.fetch_info("SLC8A1")
#fetch.fetch_info("SLC8A1")