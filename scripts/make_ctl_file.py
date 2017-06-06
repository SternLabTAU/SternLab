#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename
from PAML_utilities import write_ctl_file

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-c", "--ctl", dest="ctl", help="ctl file")
    parser.add_option("-t", "--tree", dest="tree", help="tree file")
    parser.add_option("-s", "--seq", dest="seq", help="seq file")
    parser.add_option("-o", "--out", dest="out", help="out result file")
    parser.add_option("-m", "--model", dest="model", default = "8", help="model number: 7:REV, 8:UNREST")
    parser.add_option("--fix_alpha", dest="fix_alpha", default="1", help="fix alpha")
    parser.add_option("--alpha", dest="alpha", default="0", help="alpha")

    (options, args) = parser.parse_args()
    
    ctl = options.ctl
    tree = options.tree
    seq = options.seq
    out = options.out
    model = options.model
    fix_alpha = options.fix_alpha
    alpha = options.alpha

    #check filnames
    ctl = check_filename(ctl, False)
    tree = check_filename(tree)
    seq = check_filename(seq)
    out = check_filename(out, False)
    
    ctl = write_ctl_file(ctl, seq, tree, out, model, fix_alpha, alpha)
    print("wrote ctl file to %s" % ctl)



    


if __name__ == "__main__":
    main()

