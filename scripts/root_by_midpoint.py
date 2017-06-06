#! /usr/local/python_anaconda/bin/python3.4


from optparse import OptionParser
from phylogenetic_utilities import root_at_midpoint
from file_utilities import  check_filename

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--tree", dest="tree", help="input tree file")
    parser.add_option("-o", "--out", dest="output", help="output tree file", default=None)

    (options, args) = parser.parse_args()
    tree_file = options.tree
    output = options.output

    #check filnames
    if tree_file == None:
        parser.error("you must specify a tree filename")
    tree_file = check_filename(tree_file)
    if output != None:
        output = check_filename(output, Truefile=False)

    output = root_at_midpoint(tree_file, output)
    print("rooted %s at midpoint -> %s" % (tree_file, output))




if __name__ == "__main__":
    main()

