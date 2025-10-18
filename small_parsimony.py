import argparse
import os
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
    for node in adj_list:
        for child in adj_list[node]:
            if all(c in "ACGT" for c in child):
                seq_len = len(child)
                break

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

            # Process a node only after both children have been processed (unripe)
            # Check whether both children are not ripe and process the node
            if ((not son_ripeness) and (not daughter_ripeness)):
                for nuc in T[v]["scores"]:
                    son_name = T[v]["children"][0]
                    daughter_name = T[v]["children"][1]
                    T[v]["scores"][nuc] = find_parsimony(T[son_name], T[daughter_name], nuc)
                T[v]["ripe"] = False
                
            else:
                # Add to the beginning of ripe list
                ripe_set.add(v)

        # Reset all nodes to ripe
        for v in T:
            T[v]["ripe"] = True

        # Go back down the tree
        root_nuc = min(T[root]["scores"], key=T[root]["scores"].get)
        T[root]["sequence"] += root_nuc

        set_lowest_nucs(root_nuc, root, T)

        final_score += min(T[root]["scores"].values())
        # print("Final score: ", final_score)
        # print("Root ", root)

    final_tree = format_output_dict(T)

    return final_score, final_tree

def format_output_dict(T) -> Dict[str, int]:
    # Each entry has the format ACTGATCACTA->ACTAGCTACGA:2
    output_dict = {}
    label_to_sequence: Dict[str, str] = {}

    for node, data in T.items():
        if data["sequence"]: # if it is an internal node
            label_to_sequence[node] = data["sequence"]
        else:
            label_to_sequence[node] = node # leaf node is already a sequence

    for parent, data in T.items():
        parent_seq = label_to_sequence[parent]
        for child in data["children"]:
            child_seq = label_to_sequence[child]
            dist = hammingDistance(parent_seq, child_seq)

            output_dict[f"{parent_seq}->{child_seq}"] = dist
            output_dict[f"{child_seq}->{parent_seq}"] = dist

    return output_dict

def dict_to_string(score: int, output_dict: Dict[str, int]) -> str:
    output_string = f"{score}\n"
    for k, v in output_dict.items():
        output_string += f"{k}:{v}\n"
    return output_string

def hammingDistance(pattern: str, string: str) -> int:
    d = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            d += 1
    return d

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

    score, final_dict = calculate_small_parsimony(n, edge_dict)
    output = dict_to_string(score, final_dict)

    # Create an output directory if it doesn't already exist
    if not os.path.exists("outputs"):
        os.makedirs("outputs")

    with open('outputs/parsimony_output.txt', 'w') as f:
        f.write(output)
    


if __name__ == "__main__":
    main()