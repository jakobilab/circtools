# This script is used to fetch exon sequences per circle 
# using the REST API by ENSMBLE

import sys
import requests

class sequence(object):

    def __init__(self, geneid) -> None:
        self.geneid = geneid

    def fetch_sequence(self):

        server = "https://rest.ensembl.org"
        ext = "/sequence/id/" + self.geneid + "?"
 
        r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
 
        if not r.ok:
            r.raise_for_status()
            sys.exit() 
        self.r = r
        print(r.text)

#obj = sequence("ENSCAFG00845024998")
#obj.fetch_sequence()