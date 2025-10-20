from typing import Dict, List, Tuple
import copy
import argparse

"""
Nearest Neighbors of a Tree Problem: Given an edge in a binary tree, generate the two neighbors of this tree.

    Input: An internal edge in a binary tree.
    Output: The two nearest neighbors of this tree (for the given internal edge).
"""

def nearest_neighbors(edge:List[str], tree: Dict[str, List[str]]):
    # Identify internal edge nodes
    node_A = edge[0]
    node_B = edge[1]
    print("A: ", node_A)
    print("B: ", node_B)
    print(tree)
    # Pop internal edge nodes to end of adjacency list
    tree[node_A].remove(node_B)
    tree[node_A].append(node_B)

    tree[node_B].remove(node_A)
    tree[node_B].append(node_A)

    # Create deep copies for each neighbor
    neighbor1= copy.deepcopy(tree)
    neighbor2 = copy.deepcopy(tree)

    # Identify switch node and other nodes
    # FIX: How to check pos 0 and 1 to see if they are correct
    switch_node = tree[node_A][0]
    switch1 = tree[node_B][0]
    switch2 = tree[node_B][1]

     # Forward edge
    neighbor1[switch_node] = node_B
    neighbor2[switch_node] = node_B

    # Backward edge
    neighbor1[node_B][0] = switch_node
    neighbor2[node_B][1] = switch_node

    neighbor1[switch1] = node_A
    neighbor2[switch2] = node_A

    neighbor1[node_A][0] = switch1
    neighbor2[node_A][1] = switch2

    return neighbor1, neighbor2

def parse_input(input_lines: List[str]) -> Dict[str, List[str]]:
    print("INPUT: ", input_lines[0].strip().split(" "))
    edge = input_lines[0].strip().split(" ")
    tree_dict: Dict[str, List[str]] = {}

    # Split on the -> and save everything it points to an adjacency list
    for line in input_lines[1:]:
        node, rest = line.strip().split("->")
        if node not in tree_dict:
            tree_dict[node] = [rest]
        else:
            tree_dict[node].append(rest)

    return edge, tree_dict

def print_neighbors(neighbor_tuple) -> str:
    tree1 = neighbor_tuple[0]
    tree2 = neighbor_tuple[1]

    output_string = f""

    # Formatting tree1
    for k, v in tree1.items():
        for child in v:
            output_string += f"{k}->{child}\n"

    output_string += "\n"

    # Formatting tree2
    for k, v in tree2.items():
        for child in v:
            output_string += f"{k}->{child}\n"
    return output_string

def main():
    # Parse text file
    parser = argparse.ArgumentParser(description="Process edge-weighted graph.")
    parser.add_argument("input_file", help="Path to input file.")

    args = parser.parse_args()

    with open(args.input_file, 'r') as file:
        input_data = file.readlines()

    edge, tree = parse_input(input_data)

    # Return and print two nearest neighbors
    output = print_neighbors(nearest_neighbors(edge, tree))

    print(output)
    pass

if __name__ == "__main__":
    main()