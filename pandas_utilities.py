#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
from file_utilities import check_filename


def merge_dfs(files, output):
    """
    merges several dataframe files to output file
    :param files: input dataframe files paths
    :param output:  output file path
    """
    files = [check_filename(f) for f in files]
    output = check_filename(output, Truefile=False)
    dfs = [pd.read_csv(f) for f in files]
    merged = df.append(dfs)
    df = pd.DataFrame()
    merged.to_csv(output)



if __name__ == "__main__":
    main()
