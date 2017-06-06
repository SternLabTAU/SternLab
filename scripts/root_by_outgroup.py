#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename
from phylogenetic_utilities import root_tree

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--tree", dest="tree_file", help="input tree file")
    parser.add_option("-o", "--output", dest="output_file", help="tree output file")
    parser.add_option("--id_outgroup_file", dest="id_outgroup_file", help="file that contains id for outgrouping")
    parser.add_option("--id_outgroup", dest="id_outgroup", help="id for outgrouping")

    (options, args) = parser.parse_args()
    tree_file = options.tree_file
    output = options.output_file
    outgroup_file = options.id_outgroup_file
    outgroup_id = options.id_outgroup

    # check filenames an options
    if outgroup_file == None and outgroup_id == None:
        parser.error("need to specify outgroup ID or outgroup file that contains outgroup ID")
    elif outgroup_file != None and outgroup_id != None:
        parser.error("need to specify only one - outgroup ID or outgroup file that contains outgroup ID")

    tree_file = check_filename(tree_file)
    output = check_filename(output, Truefile=False)
    if outgroup_file != None:
        outgroup_file = check_filename(outgroup_file)
        outgroup_id = open(outgroup_file, "rb").read()

    root_tree(tree_file, output, outgroup_id)


if __name__ == "__main__":
    main()