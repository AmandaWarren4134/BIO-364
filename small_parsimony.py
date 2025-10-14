import argparse
#import numpy as np
import copy
from typing import Dict, Tuple, List


def calculate_small_parsimony(n: int, adj_list: Dict[str, List[str]]) -> Tuple[int, Dict[str, List]]:
    # Convert adj_list to T
    T = build_T(adj_list)

    # Create a set of ripeness
    ripe_set = copy.deepcopy(T.keys())

    root = ""

    # Loop through each nucleotide in the sequence
    # FIX THIS LATER
    i = 0

    # Instantiate leaf nodes as 0 and infinty
    for v in T:
        if (T[v]["children"] == []):
            # Set leaf nodes as unripe and remove them from set
            T[v]["ripe"] = False
            ripe_set.remove(v)
            for nuc in T[v]["scores"]:
                if (nuc == v[i]): # DEBUG: Check if this works or if it should be T[v].tostring()[i]
                    T[v]["scores"][nuc] = 0
    
    while (ripe_set):
        # If there is one node left, label it as the root
        if (len(ripe_set) == 1):
            root = ripe_set[0]

        # If both children are unripe, then 
        if ():
            T[v]["ripe"] = False
            ripe_set.remove(v)
            for nuc in T[v]["scores"]:
                son_name = T[v]["children"][0]
                daughter_name = T[v]["children"][1]

                # Find the parsimony score for this nucleotide
                nuc_score = find_parsimony(T[son_name], T[daughter_name], nuc)
                T[v]["scores"][nuc] = nuc_score
    
    # Go back down the tree
    for v in T:
        if T[v]["children"] == []:
            pass
    
    return min(T[root]["scores"].values())  

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
    # Split on the -> and save everything it points to to an adjacency list
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

    print("N: ", n)
    print("Edge_dict: ", edge_dict)

    # tree = {
    # "4": ["CAAATCCC", "ATTGCGAC"],
    # "5": ["CTGCGCTG", "ATGGACGA"],
    # "6": ["4", "5"]
    # }
    # print(build_T(tree))

if __name__ == "__main__":
    main()