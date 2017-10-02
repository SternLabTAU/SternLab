from optparse import OptionParser
from SternLab.Candida_Positive_Selection import create_gene_table
from SternLab.Candida_Positive_Selection import create_genes_for_samples

def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-f", "--fasta", dest="FASTA_FILE_PATH", help="Reference genome fasta file")
    parser.add_option("-g", "--gene_table", dest="GENE_TABLE_PATH", help="Gene table desired location")
    parser.add_option("-p", "--patient1", dest="PATIENT1", help="Path to patient1's data")
    parser.add_option("-n", "--new_path", dest="PATH_PREFIX", help="The prefix for the gene separated fasta files")
    parser.add_option("-e", "--error_file", dest="ERRORS", help="The path to the errors file")

    (options,args) = parser.parse_args()

    FASTA_FILE_PATH = options.FASTA_FILE_PATH
    GENE_TABLE_PATH = options.GENE_TABLE_PATH
    PATIENT1 = options.PATIENT1
    PATH_PREFIX = options.PATH_PREFIX
    ERRORS = options.ERRORS

    #create_gene_table.create_gene_table(FASTA_FILE_PATH, GENE_TABLE_PATH)
    create_genes_for_samples.create_gene_files(FASTA_FILE_PATH, GENE_TABLE_PATH, PATIENT1,PATH_PREFIX, ERRORS)


if __name__ == "__main__":
    main()
