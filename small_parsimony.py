import argparse
#import numpy as np
import copy
from functools import partial
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
                son_scores = T[v]["children"][0]
                daughter_scores = T[v]["children"][1]

                # Find the parsimony score for this nucleotide
                nuc_score = find_parsimony(son_scores, daughter_scores, nuc)
                T[v]["scores"][nuc] = nuc_score
    
    return min(T[root]["scores"].values())
    

# Myesha
def find_parsimony(son_scores: Dict[str, int], daughter_scores: Dict[str, int], nuc: str) -> int:
    son = copy.deepcopy(son_scores)
    daughter = copy.deepcopy(daughter_scores)
    for n in son:
        if n != nuc:
            son[n] += 1
    for n in daughter:
        if n != nuc:
            daughter[n] += 1
    son_min = min(son.values())
    daughter_min = min(daughter.values())
    return son_min + daughter_min

def build_T(adj_list: Dict[str, List[str]]) -> Dict[str, List[Dict]]:
    T = {}
    attr_dict = {"children": [],
                 "ripe" : True,
                 "sequence" : "",
                 "scores" : {"A" : 99999, "C" : 99999, "G" : 99999, "T" : 99999}
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
def process_lines(input_lines: List[str]) -> Tuple[int, Dict[str, str]]:
    n = int(input_lines[0].strip())
    edge_dict: Dict[str, str] = {}

    for line in input_lines[1:]:
        node, rest = line.strip().split("->")
        edge_dict[node] = str()

    return n, edge_dict

def main():
    # parser = argparse.ArgumentParser(description="Process edge-weighted graph.")
    # parser.add_argument("input_file", help="Path to input file.")

    # args = parser.parse_args()

    # with open(args.input_file, 'r') as file:
    #     input_data = file.readlines()
    
    # n, edge_dict = process_lines(input_data)

    # tree = {
    # "4": ["CAAATCCC", "ATTGCGAC"],
    # "5": ["CTGCGCTG", "ATGGACGA"],
    # "6": ["4", "5"]
    # }
    # print(build_T(tree))

    scores_1 = {"A" : 1, "C" : 1, "G" : 2, "T" : 2}
    scores_2 = {"A" : 2, "C" : 1, "G" : 1, "T" : 2}
    print(find_parsimony(scores_1, scores_2, 'A'))
    print(find_parsimony(scores_1, scores_2, 'C'))
    print(find_parsimony(scores_1, scores_2, 'G'))
    print(find_parsimony(scores_1, scores_2, 'T'))

if __name__ == "__main__":
    main()