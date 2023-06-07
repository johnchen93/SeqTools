# read in a fasta file and plot the histogram of sequence lengths
import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt

def Plot(fasta, binwidth=25):
    lengths = []
    for record in SeqIO.parse(fasta, "fasta"):
        lengths.append(len(record.seq))

    plt.hist(lengths, bins=range(min(lengths), max(lengths) + binwidth, binwidth))

    # label axes
    plt.xlabel("Sequence Length")
    plt.ylabel("Count")

    plt.show()

if __name__ == "__main__":
    # use argparse for commandline interface
    parser = argparse.ArgumentParser(description="Plot the histogram of sequence lengths in a fasta file.")
    parser.add_argument("fasta", help="fasta file to plot")
    parser.add_argument("-b", "--binwidth", help="bin width for histogram", type=int, default=25)
    args = parser.parse_args()
    Plot(args.fasta)
