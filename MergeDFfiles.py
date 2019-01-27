import sys
import pandas as pd
import glob

file_list = glob.glob(sys.argv[1])
output_file = sys.argv[2]

dfs = [pd.read_csv(f, sep="\t") for f in files]
df = pd.DataFrame()
merged = df.append(dfs)

merged.to_csv(output_file, index=False, sep="\t")