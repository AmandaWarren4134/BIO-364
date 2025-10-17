import argparse
#import numpy as np
import copy
from typing import Dict, Tuple, List

# from setuptools.dist import sequence


def calculate_small_parsimony(n: int, adj_list: Dict[str, List[str]]) -> Tuple[int, Dict[str, int]]:
    # Convert adj_list to T
    T = build_T(adj_list)

    root = ""
    final_score = 0

    # Find the root
    all_children = []
    for node, child_list in adj_list.items():
        for child in child_list:
            all_children.append(child)
    
    for node in adj_list:
        if node not in all_children:
            root = node

    # Find the length of a leaf node sequence
    seq_len = 0
    for child_list in adj_list.values():
        if len(child_list[0]) > seq_len:
            seq_len = len(child_list[0])

    # Loop through each nucleotide in the sequence
    for i in range(seq_len):

        # Create a set of ripeness
        ripe_set = set(T.keys())

        # Instantiate leaf nodes as 0 and infinty
        for v in T:
            if (T[v]["children"] == []):
                # Set leaf nodes as unripe and remove them from set
                T[v]["ripe"] = False
                ripe_set.remove(v)
                for nuc in T[v]["scores"]:
                    if (nuc == v[i]):
                        T[v]["scores"][nuc] = 0
                    else:
                        T[v]["scores"][nuc] = 99999
        
        while (ripe_set):
            v = ripe_set.pop()
            # # If there is one node left, label it as the root
            # if (len(ripe_set) == 1):
            #     root = v

            son_ripeness = T[T[v]["children"][0]]["ripe"]
            daughter_ripeness = T[T[v]["children"][1]]["ripe"]

            # If both children are unripe, then 
            if ((not son_ripeness) and (not daughter_ripeness)):
                T[v]["ripe"] = False
                
                for nuc in T[v]["scores"]:
                    son_name = T[v]["children"][0]
                    daughter_name = T[v]["children"][1]

                    # Find the parsimony score for this nucleotide
                    nuc_score = find_parsimony(T[son_name], T[daughter_name], nuc)
                    #DEBUG
                    # print(f"Ripe set: ", ripe_set)
                    # print(f"Nuc parsimony score: {nuc_score}\t Current node: {v}\t Nuc: {nuc}\n")
                    T[v]["scores"][nuc] = nuc_score
            else:
                # Add to the beginning of ripe list
                ripe_set.add(v)

        # Reset all nodes to ripe
        for v in T:
            T[v]["ripe"] == True

        # Go back down the tree
        root_nuc = min(T[root]["scores"], key=T[root]["scores"].get)
        T[root]["sequence"] += root_nuc

        set_lowest_nucs(root_nuc, root, T)

        final_score += min(T[root]["scores"].values())
        # print("Final score: ", final_score)
        # print("Root ", root)

    final_tree = format_output_dict(T)

    return final_score, final_tree

# MYESHA
def format_output_dict(T) -> Dict[str, int]:
    # Each entry has the format ACTGATCACTA->ACTAGCTACGA:2

    t: Dict[str, list[str]] = {}
    output_dict = {}
    for v in T:
        if len(T[v]['sequence']) > 0: #has a sequence / not a leaf node
            children: list[str] = []
            for child in T[v]['children']:
                if len(T[child]['sequence']) > 0:
                    children.append(T[child]['sequence'])
                else:
                    children.append(child)
            t[T[v]['sequence']] = children
        else:
            t[v] = []

    for v in t:
        if len(t[v]) > 0:
            son_length = hammingDistance(v, t[v][0])
            daughter_length = hammingDistance(v, t[v][1])
            output_dict[f"{v}->{t[v][0]}"] = son_length
            output_dict[f"{t[v][0]}->{v}"] = son_length
            output_dict[f"{v}->{t[v][1]}"] = daughter_length
            output_dict[f"{t[v][1]}->{v}"] = daughter_length

    return output_dict

def hammingDistance(pattern: str, string: str) -> int:
    d = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            d += 1
    return d

# AMANDA
def dict_to_string(score: int, output_dict: Dict[str, int]) -> str:
    output_string = f"{score}\n"
    for k, v in output_dict.items():
        output_string += f"{k}:{v}\n"
    return output_string

def set_lowest_nucs(parent_nuc: str, node: str, T) -> None:

    if T[node]["children"] == []:
        return

    son_name = T[node]["children"][0]
    daughter_name = T[node]["children"][1]

    son_nuc = T[son_name]["nuc_choices"][parent_nuc]
    daughter_nuc = T[daughter_name]["nuc_choices"][parent_nuc]

    T[son_name]["sequence"] += son_nuc
    T[daughter_name]["sequence"] += daughter_nuc

    clear_nuc_choices(T[son_name]["nuc_choices"])
    clear_nuc_choices(T[daughter_name]["nuc_choices"])

    set_lowest_nucs(son_nuc, son_name, T)
    set_lowest_nucs(daughter_nuc, daughter_name, T)

def clear_nuc_choices(choices_dict: Dict[str, str]) -> None:
    for nuc in choices_dict:
        choices_dict[nuc] = ""

# Myesha
def find_parsimony(son_node, daughter_node, nuc: str) -> int:
    son = copy.deepcopy(son_node["scores"])
    daughter = copy.deepcopy(daughter_node["scores"])

    for n in son:
        if n != nuc:
            son[n] += 1

    for n in daughter:
        if n != nuc:
            daughter[n] += 1

    # Minimum score for son and daughter
    son_min = min(son.values())
    daughter_min = min(daughter.values())

    # Find best nucleotide
    son_nuc = min(son, key=son.get)
    daughter_nuc = min(daughter, key=daughter.get)

    # Set nuc choice to best nuc
    son_node["nuc_choices"][nuc] = son_nuc
    daughter_node["nuc_choices"][nuc] = daughter_nuc

    return son_min + daughter_min

def build_T(adj_list: Dict[str, List[str]]) -> Dict[str, List[Dict]]:
    T = {}
    attr_dict = {"children": [],
                 "ripe" : True,
                 "sequence" : "",
                 "scores" : {"A" : 99999, "C" : 99999, "G" : 99999, "T" : 99999}, 
                 "nuc_choices" : {"A" : "", "C" : "", "G" : "", "T" : ""}
                 }
    
    # Loop through adjacency list and create a dictionary T
    for node in adj_list:
        if node not in T.keys():
            T[node] = copy.deepcopy(attr_dict)
        for subitem in adj_list[node]:
            T[node]["children"].append(subitem)
            if subitem not in T.keys():
                T[subitem] = copy.deepcopy(attr_dict)

    return T

# AMANDA
# Take the input .txt file and parse it to grab the integer and create a dictionary
def process_lines(input_lines: List[str]) -> Tuple[int, Dict[str, List[str]]]:
    # Take the first line as n
    n = int(input_lines[0].strip())
    edge_dict: Dict[str, List[str]] = {}

    # Each line is in the format 4->ACCTGCAGCTCA
    # Split on the -> and save everything it points to an adjacency list
    for line in input_lines[1:]:
        node, rest = line.strip().split("->")
        if node not in edge_dict:
            edge_dict[node] = [rest]
        else:
            edge_dict[node].append(rest)

    return n, edge_dict

def main():
    parser = argparse.ArgumentParser(description="Process edge-weighted graph.")
    parser.add_argument("input_file", help="Path to input file.")

    args = parser.parse_args()

    with open(args.input_file, 'r') as file:
        input_data = file.readlines()

    n, edge_dict = process_lines(input_data)

    # print("N: ", n)
    # print("Edge_dict: ", edge_dict)

    score, final_dict = calculate_small_parsimony(n, edge_dict)
    print(dict_to_string(score, final_dict))

    tree = {
    "4": ["CAAATCCC", "ATTGCGAC"],
    "5": ["CTGCGCTG", "ATGGACGA"],
    "6": ["4", "5"]
    }
    #print(build_T(tree))
    # print(calculate_small_parsimony(4, tree))



if __name__ == "__main__":
    main()