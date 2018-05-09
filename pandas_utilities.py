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


def add_column_to_csv(input_df_file, new_column_name, new_column_value, output_df_file=None, sep=None):
    input_df_file = check_filename(input_df_file)
    if output_df_file != None:
        output_df_file = check_filename(output_df_file, Truefile=False)
    else:
        output_df_file = input_df_file
    if sep != None:
        df = pd.read_csv(input_df_file, sep=sep)
    else:
        df = pd.read_csv(input_df_file)
    df[new_column_name] = new_column_value
    df.to_csv(output_df_file)


if __name__ == "__main__":
    main()
