#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from pbs_runners import blast_runner, mafft_runner, njTree_runner, prank_runner, phyml_aa_runner, pear_runner
from pbs_runners import sampling_runner, umerge_runner, ufilter_runner
from file_utilities import check_filename


def main():
    """
    runs some of the pbs runners functions
    """
    parser = OptionParser("usage: %prog [options]")
    #which pbs runners funtions to run
    parser.add_option("-b", "--blast", dest="blast", action="store_true",  default=False,  help="do you want to run blast runner")
    parser.add_option("-s", "--sampler", dest="sampler", action="store_true",  default=False,  help="do you want to run sampler runner")
    parser.add_option("-p", "--prank", dest="prank", action="store_true",  default=False,  help="do you want to run prank runner")
    parser.add_option("-r", "--pear", dest="pear", action="store_true",  default=False,  help="do you want to run pear runner")
    parser.add_option("-m", "--mafft", dest="mafft", action="store_true",  default=False,  help="do you want to run mafft runner")
    parser.add_option("-n", "--njTree", dest="njTree", action="store_true",  default=False,  help="do you want to run njTree runner")
    parser.add_option("-l", "--phyml", dest="phyml", action="store_true",  default=False,  help="do you want to run phyml runner")
    parser.add_option("-u", "--umerge", dest="umerge", action="store_true",  default=False,  help="do you want to run umerge runner")
    parser.add_option("-f", "--ufilter", dest="ufilter", action="store_true",  default=False,  help="do you want to run ufilter runner")
    #options for functions
    parser.add_option("-i", "--input", dest="input_file", help="input file")
    parser.add_option("-o", "--output", dest="output_file", help="output file", default=None)
    parser.add_option("-t", "--input2", dest="input_file2", help="input file2", default=None)
    parser.add_option("-a", "--amount", dest="amount", help="amount")

    (options, args) = parser.parse_args()
    blast = options.blast
    mafft = options.mafft
    prank = options.prank
    njTree = options.njTree
    phyml = options.phyml
    pear = options.pear
    sampler = options.sampler
    umerge = options.umerge
    ufilter = options.ufilter


    input_file = options.input_file
    output_file = options.output_file
    input_file2 = options.input_file2
    amount = options.amount

    input_file = check_filename(input_file)
    if input_file2 != None:
        input_file2 = check_filename(input_file2)
    if output_file != None:
        output_file = check_filename(output_file, Truefile=False)


    if blast:
        blast_runner(input_file, outfile=output_file)
    elif mafft:
        mafft_runner(input_file, output_file)
    elif njTree:
        njTree_runner(input_file, output_file)
    elif prank:
        prank_runner(input_file, output_file)
    elif phyml:
        phyml_aa_runner(input_file)
    elif pear:
        pear_runner(input_file, input_file2, output_file)
    elif sampler:
        sampling_runner(input_file, amount, output_file)
    elif umerge:
        umerge_runner(input_file, output_file)
    elif ufilter:
        ufilter_runner(input_file, output_file)



if __name__ == "__main__":
    main()
