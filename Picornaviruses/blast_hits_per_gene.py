import pandas as pd
import xml.etree.ElementTree as ET
from optparse import OptionParser

# blast_path = '/sternadi/home/volume2/ella/Picornaviruses/samples/rhinovirus/blast analysis/blast_results.txt'
# output = '/sternadi/home/volume2/ella/Picornaviruses/samples/rhinovirus/blast analysis/'


def hits_per_gene(blast_path, gene_info, output_path, virus):
    '''
    Reads through the blast output file and sums how many hits came up for every gene
    :param blast_path: blast result file
    :param gene_info: dictionary of all gene positions
    :param output_path: common path for all result files
    :param virus: virus name - rhino, HevC, HevB, RVA-minor
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

        # confirm the hit is a of the right virus
        hit_description = hit[2].text.lower()
        if virus == "rhino" and not is_rhino(hit_description):
            continue
        elif virus == "HevC" and not is_HevC(hit_description):
            continue
        elif virus == "HevB" and not is_HevB(hit_description):
            continue
        elif virus == "RVA-minor" and not is_rva_minor(hit_description):
            continue
        elif virus == "RVA-major" and not is_rva_major(hit_description):
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


def return_HevB_genes():
    genes = {
        "VP1" : (2452,3303),
        "3D" : (5911,7296)
    }
    return genes


def return_rva_minor_genes():
    genes = {"VP1" : (2233, 3081)}
    return genes

def return_rva_major_genes():
    genes = {"VP1" : (2017, 2904)}
    return genes

def is_rva_minor(description):
    irrelevant_phrases = ["utr", "rhinovirus c", "rhinovirus b", "untranslated"]
    minor_numbers = [23, 30, 2, 49, 31, 47, 25, 62, 29, 44, "1A", "1B"]
    rhinoA_phrases = {"human rhinovirus a", "hrv-a", "hrva"}
    for irrelevant_phrase in irrelevant_phrases:
        if irrelevant_phrase in description:
            return False
    for rhinoA_phrase in rhinoA_phrases:
        if rhinoA_phrase in description:
            for minor_phrase in minor_numbers:
                if "a" + str(minor_phrase) in description:
                    return True
    return False


def is_rva_major(description):
    irrelevant_phrases = ["utr", "rhinovirus c", "rhinovirus b", "untranslated"]
    majors = ["a1", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "a15", "a16", "a18",
              "a19", "a20", "a21", "a22", "a24", "a28", "a32",
              "a33", "a34", "a36", "a38", "a39", "a40", "a41", "a43", "a45", "a46",
              "a50", "a51", "a53", "a54", "a55", "a56", "a57", "a58", "a59", "a60", "a61",
              "a63", "a64", "a65", "a66", "a67", "a68", "a71", "a73", "a74", "a75", "a76", "a77",
              "a78", "a80", "a81", "a82", "a85", "a88", "a89", "a90", "a94", "a96", "a100", "a101",
              "a102", "a103", "a104", "a105", "a106", "a107", "a108", "a109"]
    rhinoA_phrases = {"human rhinovirus a", "hrv-a", "hrva"}
    for irrelevant_phrase in irrelevant_phrases:
        if irrelevant_phrase in description:
            return False
    for rhinoA_phrase in rhinoA_phrases:
        if rhinoA_phrase in description:
            for major in majors:
                if major in description:
                    return True
    return False

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


def is_HevC(description):
    '''
    Check if the hit belongs to a relevant HevC virus
    :param description: hit description line from blast result
    :return: True \ False
    '''
    irrelevant_phrases = ["utr", "untranslated"]
    HevC_numbers = [1, 11, 13, 17, 19, 20, 21, 22, 24]
    polio = ["human poliovirus 1", "human poliovirus 2", "human poliovirus 3"]
    for irrelevant_phrase in irrelevant_phrases:
        if irrelevant_phrase in description:
            return False
    for phrase in polio:
        if phrase in description:
            return True
    Hev_starters = ["human coxsackievirus a", "cv-a"]
    for n in HevC_numbers:
        for start in Hev_starters:
            if start + str(n) + " " in description:
                return True


def is_HevB(description):
    '''
    Check if the hit belongs to a relevant HevB virus
    :param description: hit description line from blast result
    :return: True \ False
    '''
    irrelevant_phrases = ["utr", "untranslated"]
    HevB_numbers = [1,2,3,4,5,6]
    for irrelevant_phrase in irrelevant_phrases:
        if irrelevant_phrase in description:
            return False
    HevB_starters = ["human coxsackievirus b", "cv-b"]
    for n in HevB_numbers:
        for start in HevB_starters:
            if start + str(n) + " " in description:
                return True
    if "cv-a9" in description or "human coxsackievirus a9" in description:
        return True



