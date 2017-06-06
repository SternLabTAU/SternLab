#! /usr/local/python_anaconda/bin/python3.4

from os import path
from Bio import Phylo
from file_utilities import check_dirname
from file_utilities import check_filename

def root_tree(tree_file, outgroup, output):
    """
    root file using an outgroup
    :param tree_file: input tree file path
    :param outgroup: outgroup name
    :param output: output rooted file path
    :return: output file path
    """
    tree_file = check_filename(tree_file)
    output = check_filename(output, Truefile=False)
    tree = Phylo.read(tree_file, "newick")
    outgroup_clade = list(tree.find_clades(outgroup))[0]
    tree.root_with_outgroup(outgroup_clade)
    Phylo.write(tree, output_file, "newick")
    return output


def save_all_rooted_trees(tree_file, output_dir):
    """
    saves all possible rooted trees for a given tree
    :param tree_file: input tree path
    :param output_dir: output directory for trees
    :return: rooted file dictionary
    """
    # save all possible rooted tree of a given tree
    tree_file = check_filename(tree_file)
    output_dir = check_dirname(output_dir)
    basename =  path.basename(path.splitext(tree_file)[0])
    treefile_out = output_dir + "/" + basename
    tree = Phylo.read(tree_file, "newick")
    clades = list(tree.find_clades())
    out_files = []
    for clade in clades:
        if clade.name != None:
            tree.root_with_outgroup(clade)
            if tree.rooted == True:
                outfile = treefile_out +"_%s.txt" % clade.name
                Phylo.write(tree, outfile , "newick")
                out_files.append(outfile)
                out_files[clade.name]["rooted_file"] = outfile
    return out_files

def total_branch_length(tree_file):
    """
    :param tree_file: input tree path
    :return: total banch length of the tree
    """
    tree_file = check_filename(tree_file)
    tree = Phylo.read(tree_file, "newick")
    return tree.total_branch_length()

def average_branch_length(tree_file):
    """
    calculates average branch length of the tree
    :param tree_file: input tree path
    :return: average branch length
    """
    tree_file = check_filename(tree_file)
    tree = Phylo.read(tree_file, "newick")
    clades = list(tree.find_clades())
    sum = 0
    count = 0
    for clade in clades:
        if clade.branch_length == None:
            continue
        sum += clade.branch_length
        count += 1
    return float(sum)/count

def maximum_branch_length(tree_file):
    """
    finds the maximal branch length of a given tree
    :param tree_file: input tree path
    :return: maximal branch length
    """
    tree_file = check_filename(tree_file)
    tree = Phylo.read(tree_file, "newick")
    clades = list(tree.find_clades())
    maximum = 0
    for clade in clades:
        if clade.branch_length == None:
            continue
        if clade.branch_length > maximum:
            maximum = clade.branch_length
    return maximum

def branch_over_threshold(tree_file, threshold):
    """
    returns True if there is a branch over the threshold in a given tree
    :param tree_file: input tree path
    :param threshold: branch length threshold
    :return: True / False
    """
    tree_file = check_filename(tree_file)
    tree = Phylo.read(tree_file, "newick")
    clades = list(tree.find_clades())
    for clade in clades:
        if clade.branch_length >= threshold:
            return True
    return False


def root_at_midpoint(tree_file, outfile = None):
    """
    root tree at midpoint
    :param tree_file: input tree path
    :param outfile: output rooted tree path (default: None)
    :return: output rooted path
    """
    tree_file = check_filename(tree_file)
    if outfile == None:
        outfile = tree_file.split(".tree")[0] + ".rooted.tree"
    else:
        outfile = check_filename(outfile, Truefile=False)
    tree = Phylo.read(tree_file, "newick")
    rooted_tree = tree.root_at_midpoint()
    Phylo.write([tree], outfile, "newick")
    return outfile



