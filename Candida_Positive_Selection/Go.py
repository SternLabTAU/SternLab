import pandas as pd
from optparse import OptionParser

def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-a", "--annotations_association", dest="ANNOTATION", help="The annotation association file for all candida genes")
    parser.add_option("-s", "--site_summary", dest="SUMMARY", help="Positions per gene file that summarizes all genes with their positively selcted sites")
    parser.add_option("-o", "--output", dest="OUTPUT", help="Output path")
    parser.add_option("-b", "--annotations_obo", dest="OBO", help="The annotation file that maps every GO id to its description")

    (options, args) = parser.parse_args()

    ANNOTATION = options.ANNOTATION
    SUMMARY = options.SUMMARY
    OUTPUT = options.OUTPUT
    OBO = options.OBO

    df = find_ids_for_gene(SUMMARY, ANNOTATION)
    add_annotation_description(df, OBO, OUTPUT)


#df = Go.find_ids_for_gene("/sternadi/home/volume2/ella/Candida/pos_per_gene.csv", "/sternadi/home/volume2/ella/Candida/annotation_data/gene_association.cgd")
def find_ids_for_gene(gene_file, annotation_file):
    """
    Creates a DataFrame that displays every gene that was significantly positively selected and every annotation it has
    :param gene_file: Path to the file that maps every gene to the number of the position where it was positively changed
    :param annotation_file: The "gene_association" file that contains all the candida genes and their annotations
    """
    table = pd.read_csv(gene_file)
    genes = list(table["Gene"])
    gene_id_table = pd.DataFrame(columns=["Gene", "GO_ID", "Type"])
    cnt = 0
    for gene in genes:
        annotations = open(annotation_file, "r")
        for line in annotations:
            if "!" in line or "taxon:5476" not in line:
                continue
            line = line.split("\t")
            names = line[10].split("|")
            if gene in names or gene == line[2]:
                cnt += 1
                go_id = line[4]
                go_type = line[8]
                gene_id_table.loc[cnt] = [gene, go_id, go_type]
    return gene_id_table

#Go.add_annotation_description(df, "/sternadi/home/volume2/ella/Candida/annotation_data/go.obo", "/sternadi/home/volume2/ella/Candida/annotations.csv")
def add_annotation_description(df, annotations_dic, output):
    """
    This function recieves the DF from the function find_ids_for_gene, adds the description of the annoation according to the
    GO id, and writes it to an output file
    :param df: The DF to be filled
    :param annotations_dic: The obo file that mapps every GO id to a description
    :param output: The path to the desired output file
    """
    df["Anotation"] = ""
    file = open(annotations_dic, "r")
    log = open("/sternadi/home/volume2/ella/Candida/log.txt", "w")
    while True:
        line1 = file.readline()
        if not line1:
            break
        if "id: GO:" in line1:
            line2 = file.readline()
        if "id: GO:" in line1 and "name:" in line2:
            go_id = line1.split(" ")[1][:-1]
            log.write(line1 + line2 + "\n")
            name = line2[5:]
            log.write(go_id + name + "\n")
            df.loc[df['GO_ID'] == go_id, "Anotation"] = name
    file.close()
    df.to_csv(output)

# Go.extract_all_genes_GOID("/sternadi/home/volume2/ella/Candida/genes.csv","/sternadi/home/volume2/ella/Candida/annotation_data/gene_association.cgd", "/sternadi/home/volume2/ella/Candida/annotation_data/all_genes_annotations.csv")
def extract_all_genes_GOID(all_candida_genes, annotation_file, output):
    """
    :param all_candida_genes: The path to the file that contains all the candida genes with their positions
    :param annotation_file:  The "gene_association" file that contains all the candida genes and their annotations
    :param output: The file were the output will be written
    """
    genes = pd.read_csv(all_candida_genes)
    gene_list = list(genes["Gene"])
    gene_id_table = pd.DataFrame(columns=["Gene", "GO_ID", "Type"])
    gene_counter = 1
    cnt = 0
    for gene in gene_list:
        print(gene_counter)
        gene_counter += 1
        annotations = open(annotation_file, "r")
        for line in annotations:
            if "!" in line or "taxon:5476" not in line:
                continue
            line = line.split("\t")
            names = line[10].split("|")
            if gene in names or gene == line[2]:
                cnt += 1
                go_id = line[4]
                go_type = line[8]
                gene_id_table.loc[cnt] = [gene, go_id, go_type]
        annotations.close()

    gene_id_table.to_csv(output)

if __name__ == "__main__":
    main()