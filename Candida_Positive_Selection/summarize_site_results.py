import file_utilities
import pandas as pd
from optparse import OptionParser
#analyze.create_positions_table('/sternadi/home/volume2/ella/Candida/results.csv', '/sternadi/home/volume2/ella/Candida/positions_info.csv', '/sternadi/home/volume2/ella/Candida/Genes/NoRef')
#analyze.extract_info_to_dic('/sternadi/home/volume2/ella/Candida/Genes/NoRef/Ca21chr7-orf19.5502/Ca21chr7-orf19.5502-M8-Result', 'Ca21chr7', 'orf19.5502')



def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-r", "--PAML_result", dest="PAML_RESULT", help="The CSV of all genes results")
    parser.add_option("-o", "--output", dest="OUTPUT", help="Output path")
    parser.add_option("-g", "--gene_folder", dest="GENES_FOLDER",
                      help="The folder that contains all candida genes subfolders, where PAML results were saved")

    (options, args) = parser.parse_args()

    GENES_FOLDER = options.GENES_FOLDER
    OUTPUT = options.OUTPUT
    PAML_RESULT = options.PAML_RESULT

    create_positions_table(PAML_RESULT, OUTPUT, GENES_FOLDER)



def create_positions_table(PAML_results, output, path_to_PAML_output):
    """
    Iterates over the PAML result csv file and for every gene that had positive selection (Pval < 0.05),
    adds a line to a new data frame:
    Chromosome  Gene    Position    AA  Certainty
    :param PAML_results:        The PAML result CSV created by the script "process results'
    :param output:              The path to the new CSV that will contain the new DF
    :param path_to_PAML_output: The prefix of all the PAML original results files (per gene), so that
                                path_to_PAML_output + chromosome + genename + M8-Result is a valid result file.
    """
    table = pd.DataFrame(columns=["Chromosome", "Gene", "Position", "AA", "Certainty"])
    PAML_result = pd.read_csv(PAML_results)
    PAML_result = PAML_result.loc[PAML_result["Adj-Pval"] < 0.05]
    for line in PAML_result.iterrows():
        gene = line[1]["Gene"]
        chromosome = line[1]["Chromosome"]
        result_file_path = path_to_PAML_output + "/" + chromosome + "-" + gene + "/" + chromosome + "-" + gene + "-M8-Result"
        new_row = extract_info_to_dic(result_file_path, chromosome, gene)
        new_row = pd.DataFrame(new_row, columns=["Chromosome", "Gene", "Position", "AA", "Certainty"])
        table = table.append(new_row)
    table.to_csv(output, index=False)


def extract_info_to_dic(result_file, chromosome, gene):
    """
    Receives a PAML result file and extracts a dictionary of all the position that had w>1 and their additional data.
    Extracts the NEB data.
    :param result_file: The path to the PAML result file
    :return: A dictionary that contains all the positions from the PAML result of the gene given.
    """
    file_utilities.check_filename(result_file)
    positions_info = {
        "Chromosome": [],
        "Gene": [],
        "Position": [],
        "AA": [],
        "Certainty": []}
    full_text = list(open(result_file))
    line_index = 0
    while full_text[line_index] != "\tProb(w>1)  mean w\n":
        line_index += 1
    for line in full_text[line_index+1:]:
        if line == "\n":
            continue
        elif "lnL" in line:
            break
        else:
            line = line.split(" ")
            start = 0
            while line[start] == "":
                start += 1
            positions_info["Chromosome"] += [chromosome]
            positions_info["Gene"] += [gene]
            positions_info["Position"] += [line[start]]
            positions_info["AA"] += [line[start + 1]]
            positions_info["Certainty"] += [line[start + 7]]
    return positions_info


if __name__ == "__main__":
    main()