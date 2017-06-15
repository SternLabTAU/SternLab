#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
from file_utilities import check_filename


def merge_dfs(files, output, separator=","):
    """
    merges several dataframe files to output file
    :param files: input dataframe files paths
    :param output:  output file path
    :param sep: separator for csv (default: ,)
    """
    files = [check_filename(f) for f in files]
    output = check_filename(output, Truefile=False)
    dfs = [pd.read_csv(f, sep=separator) for f in files]
    df = pd.DataFrame()
    merged = df.append(dfs)
    merged.to_csv(output, index=False)



if __name__ == "__main__":
    main()
