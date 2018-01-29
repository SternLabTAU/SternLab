import pandas as pd
import xml.etree.ElementTree as ET
from optparse import OptionParser

import Picornaviruses.blast_hits_per_gene as blast
import Picornaviruses.check_date_availability as date
import Picornaviruses.download_seq_by_id as download
import pbs_runners

def main():

    ### sum up blast results - create a table that sums up how many hits came up for every gene

    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-b", "--blast", dest="blast_path", help="blast result file")
    parser.add_option("-o", "--output", dest="output", help="Output path")
    parser.add_option("-v", "--virus", dest="virus", help="virus type")

    (options, args) = parser.parse_args()

    blast_path = options.blast_path
    output = options.output
    virus = options.virus

    if virus == "rhino":
        gene_info = blast.return_rhino_genes()
    elif virus == "HevC":
        gene_info = blast.return_HevC_genes()
    blast.hits_per_gene(blast_path, gene_info, output)
    print("done")

#    ### check date availability for records in the blast sum up file
#
#    parser = OptionParser("usage: %prog[options]")
#    parser.add_option("-g", "--gene_hits", dest="gene_hits", help="gene hits file")
#    parser.add_option("-o", "--output", dest="output", default="same as input", help="Output path. When left empty rewrited input")
#
#    (options, args) = parser.parse_args()
#
#    gene_hits = options.gene_hits
#    if options.output == "same as input":
#        output = gene_hits
#    else:
#        output = options.output
#
#    date.add_collection_date(gene_hits, output)
#
#    ### download sequences by id and create united fasta file
#
#    parser = OptionParser("usage: %prog[options]")
#    parser.add_option("-g", "--gene_hits", dest="gene_hits", help="gene hits file")
#    parser.add_option("-o", "--output", dest="output", help="common output path. Will be completed by the id name")
#    parser.add_option("-u", "--united", dest="united", help="path for united file")
#
#    (options, args) = parser.parse_args()
#    gene_hits = options.gene_hits
#    output = options.output
#    united = options.united
#
#    download.download_seqs(gene_hits, output)
#    download.unite_seqs_to_file(output, united)
#
#    ### run prank on united file
#
#    pbs_runners.prank_codon_runner(united)
#
#    ### build tree for the alignment using phyml
#
#    parser = OptionParser("usage: %prog[options]")
#    parser.add_option("-a", "--alignment", dest="alignment", help="alignment file")
#
#    (options, args) = parser.parse_args()
#    alignment = options.alignment
#
#    pbs_runners.phyml_runner(alignment)