
"""
@Author: daniellem1

"""

'''preprocess .freqs files in order to get for each genome position it's num of reads'''



import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from cirseq_utilities import *


def main():
#RV
    # freqs_file_path = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/RVB14p2.freqs"
    # out_dir = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/plots/"
    # coverage_graph(freqs_file_path, out_dir)
#PV
    freqs_file_path = "/Volumes/STERNADILABTEMP$/volume1/okushnir/Cirseq/PV/Mahoney/P3/20170907_q23r2_blastn/PV-p3.1036617.freqs"
    out_dir = "/Volumes/STERNADILABTEMP$/volume1/okushnir/Cirseq/PV/Mahoney/P3/20170907_q23r2_blastn/plots"
    coverage_graph(freqs_file_path, out_dir)


def coverage_graph(freqs, out_dir):
    data = parse_reads(freqs)
    pos = data[0]
    reads = data[1]
    graph = plt.plot(pos, reads, color="DarkOrchid")

    plt.xlabel("Position In The Genome[bp]", fontsize=20)
    plt.ylabel("Number Of Reads", fontsize=20)
    plt.title("Coverage", fontsize=30)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim(0, (max(pos)+10))
    plt.ylim(1, (max(reads)+1000000))
    sns.set_style("darkgrid")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(out_dir + "/coverage.png", dpi=680)
    plt.close('all')
    return graph

if __name__ == "__main__":
    main()







