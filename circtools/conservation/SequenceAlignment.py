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

    def __init__(self, input_fasta, source_species, circle_name) -> None:
        self.fasta_file = input_fasta
        self.source_species = source_species
        temp_name = circle_name.split("_")
        self.circle_name = temp_name[0] + "(" + temp_name[1] + ":" + temp_name[2] + "-" + temp_name[3] + ")"

    def alignment_to_distance_matrix(self):
        # convert the output from mafft into a distance matrix
        #print(self.out_fasta)
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
                #print(leaf.name, string)
                return(string)

        self.run_mafft()

        tree_file = self.fasta_file + ".tree"
        out_png = self.fasta_file.replace(".fasta", ".png") 
        
        tree = Phylo.read(tree_file, "newick")               # this is an output file from mafft_cline() function with --treeout option
        fig = plt.figure(figsize=(11, 10), dpi=100)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False, label_func=get_label)
        fig.text(0.50, 0.02, 'Genome Versions: Human(hs)-hg38, Mouse(mm)-mm39, Pig(ss)-SusScr11, Rat(rn)-Rn7 and Dog(cl)-CanFam6', horizontalalignment='center', wrap=True)
        axes.get_yaxis().set_visible(False)
        axes.spines['top'].set_visible(False)
        axes.spines['right'].set_visible(False)
        axes.spines['left'].set_visible(False)
        plt.title("Sequence conservation tree for circle: " + self.circle_name)
        plt.xticks(fontsize=12)
        plt.show()
        plt.savefig(out_png)

    def pairwise_alignment(self):
        # perform pairwise alignment if the flag is on

        # read in the fasta file and store the length of sequence for normalising
        fasta_sequences = SeqIO.parse(open(self.fasta_file),'fasta')
        for fasta in fasta_sequences:
            if self.source_species in fasta.id:
                length = int(len(fasta.seq))

        # perform combinations and their alignment scores
        combi = combinations(SeqIO.parse(self.fasta_file , "fasta"), 2)
        aligner = Align.PairwiseAligner()
        plot_dict = {}
        for pair in combi:
            species_1 = pair[0].id.split("(")[0]
            species_2 = pair[1].id.split("(")[0]
            if ((self.source_species == species_1) or (self.source_species == species_2)):
                alignments = aligner.align(pair[0].seq, pair[1].seq)
                #print(species_1, species_2, pair[0].id, pair[1].id, alignments.score)
                plot_dict[species_1+"_"+species_2] = float(alignments.score)
        
        # normalise the alignment scores by length of sequence of source species circle
        plot_dict = {k: v / length for k, v in plot_dict.items()}
        #print(plot_dict)

        # plot as a bar plot
        species = list(plot_dict.keys())
        scores = list(plot_dict.values())
        out_bar = self.fasta_file.replace(".fasta", "_pairwise.png") 
        fig = plt.figure(figsize = (10, 10))
        plt.bar(species, scores, color ='blue', width = 0.4)
        plt.xlabel("Pairwise alignments", fontsize = 14)
        plt.ylabel("Alignment scores", fontsize = 14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title("Pairwise alignment scores for circle: " +  self.circle_name)
        plt.show()
        plt.savefig(out_bar)
'''
if __name__ == "__main__":
    obj = Alignment("/scratch/circtools2/circtools/sample_data/temp/alignment_UXS1_2_106145190_106166083_-.fasta", "hs")
    obj.draw_phylo_tree()
    obj.pairwise_alignment()
'''
