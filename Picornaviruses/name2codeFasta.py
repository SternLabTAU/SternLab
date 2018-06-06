from Bio import SeqIO
import pandas as pd


def name2code(fasta_file, coded_output, code_file):
    i = 1
    code = {}
    o = open(coded_output, "w")
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id_code = "seq" + str(i) + "$"
            code[record.id] = id_code
            i += 1
            o.write(">")
            o.write(id_code + "\n")
            o.write(str(record.seq) + "\n")

    df = pd.DataFrame.from_dict(code, orient='index')
    df.to_csv(code_file)

    handle.close()


