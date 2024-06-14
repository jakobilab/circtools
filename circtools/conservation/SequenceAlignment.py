from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib.pyplot as plt

class Alignment(object):

    def __init__(self, input_fasta) -> None:
        self.fasta_file = input_fasta


    def run_mafft(self):
        # on the input fasta sequence file, run mafft commandline
        # that stores tree and alignedment in clustalw format

        mafft_cline = MafftCommandline(input=self.fasta_file, clustalout=True, treeout=True)
        stdout, stderr = mafft_cline()

        #with open("aligned.fasta", "w") as handle:
        #handle.write(stdout)
        #align = AlignIO.read("aligned.fasta", "fasta")

    def alignment_to_distance_matrix(self):
        # convert the output tree into distance matrix

        aln = AlignIO.read('aligned.fasta', 'phylip')
        print(aln)

        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        print(dm)

        
    def draw_phylo_tree(self):
        # visulalisation of alignment results in form of phylogenetic tree    

        tree = Phylo.read("msa.fasta.tree", "newick")               # this is an output file from mafft_cline() function with --treeout option
        fig = plt.figure(figsize=(10, 10), dpi=100)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False)
        plt.show()
        plt.savefig('test.png')