import pandas as pd


def code2name(tree, code, uncoded_output):
    tree_file = open(tree, "r")
    code_file = open(code, "r")
    tree_str = tree_file.read()
    for line in code_file:
        original = str(line.split(" ")[0])
        coded = str(line.split(" ")[1])[:-1]
        tree_str = tree_str.replace(coded, original)
    o = open(uncoded_output, "w")
    o.write(tree_str)



