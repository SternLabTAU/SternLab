import pandas as pd
import xml.etree.ElementTree as ET
from optparse import OptionParser

# blast_path = '/sternadi/home/volume2/ella/Picornaviruses/samples/rhinovirus/blast analysis/blast_results.txt'
# output = '/sternadi/home/volume2/ella/Picornaviruses/samples/rhinovirus/blast analysis/'


def hits_per_gene(blast_path, gene_info, output_path):
    '''
    Reads through the blast output file and sums how many hits came up for every gene
    :param blast_path: blast result file
    :param gene_info: dictionary of all gene positions
    :param output_path: common path for all result files
    '''

    # create data sets for all genes and counters for how many hits every gene had
    res = {}
    for gene in gene_info.keys():
        data = pd.DataFrame(columns=["id", "description", "score", "from", "to", "match"])
        res[gene] = [data, 0]

    # start reading the XML file
    xml = ET.parse(blast_path)
    root = xml.getroot()
    for hit in root[8][0][4]:

        # confirm the hit is a rhinoB hit
        hit_description = hit[2].text.lower()
        if not is_rhino(hit_description):
            continue

        hit_start = int(hit[5][0][4].text)
        hit_end = int(hit[5][0][5].text)
        hit_id = hit[1].text
        hit_score = hit[5][0][2].text

        for gene in gene_info.keys():
            gene_start = gene_info[gene][0]
            gene_end = gene_info[gene][1]
            gene_len =  gene_end - gene_start
            if hit_start <= gene_start and hit_end >= gene_end:
                match = "full"
            elif (hit_start <= gene_start and hit_end >= gene_start + gene_len/2) or \
                    (hit_end >= gene_end and hit_start <= gene_start + gene_len/2) or \
                    (hit_start >= gene_start and hit_end <= gene_end and hit_end-hit_start > gene_len/2):
                match = "half"
            elif gene.lower() in hit_description:
                match = "by name"
            else:
                continue
            new_line = pd.DataFrame({"id": [hit_id],
                                     "description": [hit_description],
                                     "score": [hit_score],
                                     "from": [hit_start],
                                     "to": [hit_end],
                                     "match": [match]},
                                      columns = ["id", "description", "score", "from", "to", "match"])
            res[gene][0] = res[gene][0].append(new_line)
            res[gene][1] += 1


    # write results to csv
    counts = pd.DataFrame(columns=["gene", "hits"])
    for gene in res.keys():
        res[gene][0].to_csv(output_path + gene + ".csv")
        new_line = pd.DataFrame({"gene": [gene],
                                 "hits": [res[gene][1]]})
        counts = counts.append(new_line)
    counts.to_csv(output_path + "hits_per_gene.csv")


def return_rhino_genes():
    '''
    Returns a dictionary of all rhino genes with their start and end locations
    :return: genes dictionary
    '''
    genes = {
        "VP4" : (629, 835),
        "VP2" : (836, 1621),
        "VP3" : (1622, 2329),
        "VP1": (2330, 3196),
        "P2-A": (3197, 3634),
        "P2-B": (3635, 3925),
        "P2-C": (3926, 4915),
        "P3-A": (4916, 5170),
        "VPg": (5171, 5239),
        "3C": (5240, 5785),
        "3D": (5786, 7165),
    }
    return genes


def return_HevC_genes():
    genes = {
        "VP1" : (2485, 3402),
        "3D" : (6004, 7386)
    }
    return genes


def is_rhino(description):
    '''
    Checks if a result is a rhino B result according to the hit's description line
    :param description: hit description line from blast result
    :return: True or False
    '''
    irrelevant_phrases = ["utr", "rhinovirus c", "rhinovirus a", "untranslated"]
    rhino_phrases = ["human rhinovirus b", "hrv-b", "hrvb"]
    rhino_B_numbers = [52,104,17,70,91,69,48,100,101,102,35,83,92,79,14,72,3,6,103,37,86,26,4,5,42,99,27,93,977,84]

    for irrelevant_phrase in irrelevant_phrases:
        if irrelevant_phrase in description:
            return False
    for rhino_phrase in rhino_phrases:
        if rhino_phrase in description:
            return True
    for num in rhino_B_numbers:
        name = "rhinovirus " + str(num) + " "
        name_2 = "rhinovirus " + str(num) + ","
        name_3 = "rhinovirus type " + str(num) + ","
        name_4 = "rhinovirus type " + str(num) + " "
        if name in description or name_2 in description or name_3 in description or name_4 in description:
            return True
    return False

if __name__ == "__main__":
    main()