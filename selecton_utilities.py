#! /usr/local/python_anaconda/bin/python3.4

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re

def extract_selecton_final_params_single_file(file):
    data = open(file, "r").readlines()
    parameters = {}
    for line in data:
        if "Selecton Optimized Parameters" in line or "wYangModel Parameters:" in line or "wrtModel Parameters:" in line:
            continue
        parameter = line.split(":")[0].replace(" ", "_")
        value = line.split(":")[1].strip()
        parameters[parameter] = value


    return parameters



def extract_selecton_final_params(files, output_file=None):
    df = pd.DataFrame()

    for f in files:
        filename = f.split("/")[-1].split("_output")[0]
        model = f.split("/")[-2]
        parameters = extract_selecton_final_params_single_file(f)
        parameters["model"] = model
        parameters["protein"] = filename
        df = df.append(parameters, ignore_index=True)

    print(df)
    if output_file!=None:
        output_file = check_filename(output_file, Truefile=False)
        df.to_csv(output_file)
    return df