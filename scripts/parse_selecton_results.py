#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename, check_dirname
import glob

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--directory", dest="dir", help="input directory with kaks gaps files")
    parser.add_option("-o", "--output", dest="output", help="output file name")
    parser.add_option("-v", "--virus", dest="virus", help="virus - tells the script how to parse filenames. options: influenza, tilv")
    (options, args) = parser.parse_args()


    dir = options.dir
    output = options.output
    virus = options.virus

    dir = check_dirname(dir)
    output = check_filename(output, Truefile=False)
    files = glob.glob("%s/*kaks*gaps" % dir)
    if len(files) == 0:
        raise Exception("No files in %s" % dir)

    if virus == "influenza" or virus == "tilv": #add more virus name if you add anything
        output_text = "POS\tAMINO\tKaKs\tconfidence_interval\tvirus\tprotein\n"
    else:
        output_text = "POS\tAMINO\tKaKs\tconfidence_interval\n"

    for f in files:
        print(f)
        if virus == "influenza":
            virus_name = "Influenza " + f.split("/")[-1].split("inf")[1].split("_")[0]
            protein = f.split("_")[2]
        elif virus == "tilv":
            virus_name = "TiLV"
            protein = "Segment " + f.split("/")[-1].split("_")[0].split("seg")[-1]
        #if adding more virus type this is what you have to do:
        #elif virus == "XXXX":
            #virus_name = XXX
            #protein = XXX
        else:
            virus_name = None
            protein = None
        output_text += kaks_file_to_txt_delimeted_results(f, virus_name, protein)

    output = open(output, "w")
    output.write(output_text)
    output.close()


def kaks_file_to_txt_delimeted_results(f, virus=None, protein=None):
    kaks = open(f, "r").readlines()
    out = ""
    for line in kaks[6:]:
        if virus == None and protein == None:
            out += "\t".join(line.split("\t")[:4]) + "\n"
        else:
            out += "\t".join(line.split("\t")[:4]) + "\t%s\t%s\n" % (virus, protein)

    return out


if __name__ == "__main__":
    main()

