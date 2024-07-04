from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib.pyplot as plt
from Bio import Phylo
from itertools import combinations
from Bio import SeqIO,  Align
from Bio.Seq import Seq
from Bio.Align import *

class Alignment(object):

    def __init__(self, input_fasta, source_species) -> None:
        self.fasta_file = input_fasta
        self.source_species = source_species

    def alignment_to_distance_matrix(self):
        # convert the output from mafft into a distance matrix
        print(self.out_fasta)
        aln = AlignIO.read(self.out_fasta, 'clustal')
        #print(aln)

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

        # defining function for label drawing in phylo plot
        def get_label(leaf):
            if leaf.name == None:
                return(leaf.name )
            else:
                temp = leaf.name.split("_")
                string = temp[1] + "(" + temp[2] + ":" + temp[3] + ")"
                print(leaf.name, string)
                return(string)

        self.run_mafft()

        tree_file = self.fasta_file + ".tree"
        out_png = self.fasta_file.replace(".fasta", ".png") 
        
        tree = Phylo.read(tree_file, "newick")               # this is an output file from mafft_cline() function with --treeout option
        fig = plt.figure(figsize=(15, 10), dpi=100)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False, label_func=get_label)
        fig.text(0.50, 0.02, 'Genome Versions: Human-hg38, Mouse-mm39, Pig-SusScr11, Rat-Rn7 and Dog-CanFam6', horizontalalignment='center', wrap=True)
        axes.get_yaxis().set_visible(False)
        plt.show()
        plt.savefig(out_png)

    def pairwise_alignment(self):
        # perform pairwise alignment if the flag is on
        
        combi = combinations(SeqIO.parse(self.fasta_file , "fasta"), 2)
        aligner = Align.PairwiseAligner()
        plot_dict = {}
        for pair in combi:
            species_1 = pair[0].id.split("(")[0]
            species_2 = pair[1].id.split("(")[0]
            if ((self.source_species == species_1) or (self.source_species == species_2)):
                alignments = aligner.align(pair[0].seq, pair[1].seq)
                print(species_1, species_2, pair[0].id, pair[1].id, alignments.score)
                plot_dict[species_1+"_"+species_2] = float(alignments.score)
        
        print(plot_dict)
        # plot as a bar plot
        species = list(plot_dict.keys())
        scores = list(plot_dict.values())
        out_bar = self.fasta_file.replace(".fasta", "_pairwise.png") 
        fig = plt.figure(figsize = (10, 10))
        plt.bar(species, scores, color ='blue', width = 0.4)
        plt.xlabel("Pairwise alignments")
        plt.ylabel("Alignment scores")
        plt.show()
        plt.savefig(out_bar)

if __name__ == "__main__":
    obj = Alignment("/scratch/circtools2/circtools/sample_data/temp/alignment_UXS1_2_106145190_106166083_-.fasta", "hs")
    obj.draw_phylo_tree()
    obj.pairwise_alignment()