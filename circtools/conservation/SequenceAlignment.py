from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib.pyplot as plt
from Bio import Phylo

class Alignment(object):

    def __init__(self, input_fasta) -> None:
        self.fasta_file = input_fasta

    def alignment_to_distance_matrix(self):
        # convert the output from mafft into a distance matrix
        print(self.out_fasta)
        aln = AlignIO.read(self.out_fasta, 'clustal')
        print(aln)

        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        print(dm)

    def run_mafft(self):
        # on the input fasta sequence file, run mafft commandline
        # that stores tree and alignedment in clustalw format

        self.out_fasta = self.fasta_file + ".aligned"

        mafft_cline = MafftCommandline(input=self.fasta_file, clustalout=True, treeout=True)
        stdout, stderr = mafft_cline()

        with open(self.out_fasta, "w") as handle:
            handle.write(stdout)

        self.alignment_to_distance_matrix()

    def draw_phylo_tree(self):
        # visulalisation of alignment results in form of phylogenetic tree    

        self.run_mafft()

        tree_file = self.fasta_file + ".tree"
        out_png = self.fasta_file.replace(".fasta", ".png") 
        
        tree = Phylo.read(tree_file, "newick")               # this is an output file from mafft_cline() function with --treeout option
        fig = plt.figure(figsize=(10, 10), dpi=100)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False)
        plt.show()
        plt.savefig(out_png)

if __name__ == "__main__":
    obj = Alignment("/scratch/circtools2/circtools/sample_data/temp/alignment_UXS1_2_106145190_106166083_-.fasta")
    obj.draw_phylo_tree()