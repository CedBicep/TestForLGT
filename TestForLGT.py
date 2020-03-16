#!/usr/bin/python

# Read the bootstrapped consensus tree and ask whether the candidate LGT has phylogenetic support
# by checking various tree-based criteria for LGTs
# Arguments: python TestForLGT.py <treefile>

import csv
import os
import re
import sys
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace

target_sequence_tag = 'xxx'     # To identify candidate LGT
group_assignments = {}          # Setup group assignments
eukaryote_seqs = []             # setup a list of eukaryotic sequences in the tree
bacteria_seqs = []              # setup a list of bacterial sequences in the tree
archaea_seqs = []               # setup a list of archaeal sequences in the tree
other_seqs = []                 # setup a list of other sequences in the tree
name2tax = {}
target_leaf = ''                # setup a target leaf object
euk_supergroups = []


# Update the list of Eukaryotic supergroup : 'Viridiplantae','Oxymonadida','Alveolata'
def getAllEukSuperGroup(file):
    global euk_supergroups
    inh = open(file)
    for line in inh:
        euk_supergroups.append(line.rstrip())
    inh.close()


# Return a tree object or quit if empty
def checkInputTree(file):
    # Lets check tree string for sanity first
    inh = open(file)
    treestring = inh.readline()
    treestr = treestring.replace(';', '')
    treestr = treestr + ";"
    inh.close()
    if len(treestr) == 0:
        print(sys.argv[1] + "\tEmpty tree")
        quit()
    return Tree(treestr)


# Assign a group to each sequence of the tree
def getAnnotationForTrees(file):
    global group_assignments
    inh = open(file)
    for line in inh:
        fields = re.split('\s+', line.rstrip())
        if len(fields) >= 2:
            group_assignments[fields[0]] = fields[1]  # key = sequence ID, value = group assignment (e.g. Viridiplantae)


# Return a fully annotated tree by adding domains to each node
def annotateTree(tree):
    global name2tax
    global target_leaf
    for node in tree:
        node.add_features(domain="Other")
    for leaf in tree:
        name2tax[leaf.name] = str(str(group_assignments[leaf.name]) + "_" + "_".join(leaf.name.split('_')[:2]))
        if re.search(target_sequence_tag, leaf.name):
            leaf.add_features(domain="Eukaryote")
            eukaryote_seqs.append(leaf.name)
            target_leaf = leaf
        elif leaf.name in group_assignments:
            if group_assignments[leaf.name] in euk_supergroups:
                eukaryote_seqs.append(leaf.name)
                leaf.add_features(domain="Eukaryote")
            else:
                leaf.add_features(domain="Other")
                if group_assignments[leaf.name] in ["Bacteria"]:
                    bacteria_seqs.append(leaf.name)
                elif group_assignments[leaf.name] in ["Archaea"]:
                    archaea_seqs.append(leaf.name)
                else:
                    other_seqs.append(leaf.name)
        else:
            leaf.add_features(domain="Other")
    return tree


# Root the tree based on the biggest prokaryotic clade
def rootTreeOnProkaryotes(tree):
    # root the tree on a clade (the biggest?) of bacteria, to avoid ridiculous problems with arbitrary roots on trees
    biggest_other_node = 0
    for node in tree.get_monophyletic(values=['Other'], target_attr="domain"):
        if len(node) > biggest_other_node:
            biggest_other_node = len(node)
            tree.set_outgroup(node)
    return tree


# Test various phylogenetic criteria for LGT and write the report
def testCriteriaForLGT(tree):
    if len(eukaryote_seqs) == 1:  # this is, I guess, an LGT candidate
        Output_writer.writerow([sys.argv[1], "Singleton", "1", "N/A", "N/A", "N/A", "1"])
    else:
        try:
            answer = tree.check_monophyly(values=eukaryote_seqs, target_attr="name")
            if answer[0]:
                ca = tree.get_common_ancestor(eukaryote_seqs)
                target_group_sgs = {}
                for leaf in ca:
                    if leaf.name in group_assignments:
                        leaf_supergroup = group_assignments[leaf.name]
                        if leaf_supergroup in euk_supergroups:
                            target_group_sgs[leaf_supergroup] = 1
                    else:
                        print("Warning: a sequence in this tree do not have a supergroup assignment: " + str(leaf.name))
                num_sgs = len(target_group_sgs.keys())
                Output_writer.writerow([sys.argv[1], "Euks monophyletic", str(len(eukaryote_seqs)), str(ca.support),
                                          "N/A", "N/A", str(num_sgs)])
            elif answer[0] is False:
                mono_groups = []
                target_group = ''
                for node in tree.get_monophyletic(values=['Eukaryote'], target_attr="domain"):
                    for leaf in node:
                        if leaf.name == target_leaf.name:
                            target_group = node
                    else:
                        mono_groups.append(node)
                size_target_group = len(target_group)
                # get distance
                shortest_distance = 999999999999999.0
                # closest_other_group = ''                  # If need to do something with it
                for subtree in mono_groups:
                    curr_distance = tree.get_distance(target_group, subtree, topology_only=True)
                    if curr_distance < shortest_distance:
                        shortest_distance = curr_distance
                        # closest_other_group = subtree     # If need to do something with it
                # find out what supergroups of eukaryotes are represented in the target group
                target_group_sgs = {}
                tg_names = []
                for leaf in target_group:
                    tg_names.append(leaf.name)
                    if leaf.name in group_assignments:
                        leaf_supergroup = group_assignments[leaf.name]
                        if leaf_supergroup in euk_supergroups:
                            target_group_sgs[leaf_supergroup] = 1
                    else:
                        print("Warning: a sequence in this tree do not have a supergroup assignment: " + str(leaf.name))
                num_sgs = len(target_group_sgs.keys())
                c_a = tree.get_common_ancestor(tg_names)
                Output_writer.writerow([sys.argv[1], "Euks not monophyletic", str(len(eukaryote_seqs)),
                                          str(c_a.support), str(size_target_group), str(shortest_distance),
                                          str(num_sgs)])
            else:
                print(sys.argv[1] + "\t" + answer[0])
        except IndexError:
            raise


# Generate a colored PDF of the tree
def generateColoredPDF(tree):
    out_tree = os.path.splitext(os.path.basename(sys.argv[1]))[0] + ".pdf"
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_support = True
    ts.show_branch_length = False
    for n in tree.traverse():
        if n.name in eukaryote_seqs:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "red"
            nstyle["size"] = 15
            n.set_style(nstyle)
            n.name = name2tax.get(n.name)
        elif n.name in bacteria_seqs:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "blue"
            nstyle["size"] = 13
            n.set_style(nstyle)
            n.name = name2tax.get(n.name)
        elif n.name in archaea_seqs:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "green"
            nstyle["size"] = 13
            n.set_style(nstyle)
            n.name = name2tax.get(n.name)
        elif n.name in other_seqs:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "grey"
            nstyle["size"] = 13
            n.set_style(nstyle)
            n.name = name2tax.get(n.name)
    tree.render(out_tree, tree_style=ts)


if __name__ == "__main__":
    InputTree = checkInputTree(sys.argv[1])
    if os.path.exists('Output_file.csv'):
        append_write = 'a' # append if already exists
        Output_file = open('Output_file.csv', mode=append_write)
        Output_writer = csv.writer(Output_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    else:
        append_write = 'w' # make a new file if not
        Output_file = open('Output_file.csv', mode=append_write)
        Output_writer = csv.writer(Output_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        Output_writer.writerow(["Tree", "tResult", "EuksInTree", "SupportEukMonophyly", "EuksInTargetGroup",
                              "DistanceToClosestEukClade", "SupergroupsInTargetGroup"])
    getAnnotationForTrees("Annotation_file_for_trees.txt")
    getAllEukSuperGroup("List_that_matters.txt")
    AnnotatedTree = annotateTree(InputTree)
    RootedTree = rootTreeOnProkaryotes(AnnotatedTree)
    testCriteriaForLGT(RootedTree)
    generateColoredPDF(RootedTree)
    Output_file.close()

