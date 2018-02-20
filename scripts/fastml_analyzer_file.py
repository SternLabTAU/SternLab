#! /usr/local/python_anaconda/bin/python3.4


import os
from optparse import OptionParser

import pandas as pd

from context_analysis.fastml_analyzer import run_on_basename


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-b", "--basename", dest="basename", help="basename ")
    (options, args) = parser.parse_args()
    basename = options.basename
    
    out_file = basename + ".fastml_analysis_output.csv"
    if os.path.exists(out_file):
        print "file already exists"
        
    
    
    df = pd.DataFrame(columns = ["family", "group", "mutation", "mutation_count_in_context", "context_count_overall",
                                 "mutation_count_overall", "close_mutation"])
    df = run_on_basename(basename, df)
 
    df.to_csv(out_file)    







if __name__ == "__main__":
    main()
