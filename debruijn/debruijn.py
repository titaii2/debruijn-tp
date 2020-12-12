#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

import argparse
import os
import sys
import random
import statistics
import re
from random import randint
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.dag import descendants, ancestors
random.seed(9001)

__author__ = "Alexis Gouthey"
__copyright__ = "Université de Cergy Pontoise (CyTech)"
__credits__ = ["Alexis Gouthey"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Alexis Gouthey"
__email__ = "goutheyale@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """
    Check if path is an existing file.

    Parameters:
        path: (String) Path to the file

    Returns: The path of the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """
    Retrieve the arguments of the program.

    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest = 'fastq_file', type = isfile,
                        required = True, help="Fastq file")
    parser.add_argument('-k', dest = 'kmer_size', type = int,
                        default = 21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest = 'output_file', type = str,
                        default = os.curdir + os.sep + "contigs.fasta",
                        help = "Output contigs in fasta file")

    return parser.parse_args()


def read_fastq(fastq_file):
    """
    Read the file and return every line that contain nucleic sequences.

    Parameters:
        path: (String) Path to the fastq_file

    Returns: A generator of nucleic sequences
    """

    # Open the file in "read" mode
    file = open(fastq_file, "r")

    for line in file:
        # Match the line with regular expression containing nucleic sequences characters
        if re.match(r'[ATCG]+', line):
            # Remove the "\n" at the end of the line
            yield line.rstrip("\n")

    file.close()

def cut_kmer(read, kmer_size):
    """
    Cut a nucleic sequence with a kmer_size and return every cut sequences.

    Parameters:
        read: String of a nucleic sequence
        kmer_size: Size (Int) of returned kmer

    Returns: A generator of nucleic sequences with a size of kmer_size
    """

    for i in range(len(read) - kmer_size + 1):
        yield read[i : i + kmer_size]

def build_kmer_dict(fastq_file, kmer_size):
    """
    Build and return a dictionary containing every nucleic sequences
    with their relative number of occurence in the fastq_file.

    Parameters:
        fastq_file: Generator of nucleic sequences
        kmer_size: Size (Int) of returned kmer

    Returns: A generator of nucleic sequences and occurences with a size of kmer_size
    """

    dict = {}

    for line in read_fastq(fastq_file):
        for kmer in cut_kmer(line, kmer_size):
            if kmer in dict:
                dict[kmer] += 1
            else:
                dict[kmer] = 1
    return dict

def build_graph(kmer_dict):
    """
    Build and return a DiGraph from a kmer dictionary.

    Parameters:
        kmer_dict: Dictionary with kmer and their number of occurence

    Returns: A DiGraph
    """

    graph = nx.DiGraph()

    for kmer in kmer_dict:
        # Add edge from the dictionary
        graph.add_edge(kmer[:-1], kmer[1:], weight = kmer_dict[kmer])

    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
    Remove paths from a graph.

    Parameters:
        graph: A DiGraph
        path_list: A list of paths to delete
        delete_entry_node: (Boolean) Delete entry nodes if True
        delete_sink_node : (Boolean) Delete sink nodes if True

    Returns: A DiGraph without paths contained in path_list
    """

    for path in path_list:
        if delete_entry_node and graph.has_node(path[0]):
            graph.remove_node(path[0])

        for i in range(1, len(path) - 1):
            if graph.has_node(path[i]):
                graph.remove_node(path[i])

        if delete_sink_node and graph.has_node(path[-1]):
            graph.remove_node(path[-1])

    return graph

def std(data):
    """
    Compute the standard deviation of a list.

    Parameters:
        data: list of values (int, float)

    Returns: The standard deviation of the list "data"
    """

    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node = False, delete_sink_node = False):
    """
    Select the best path out of the path_list list and removes all others paths from the graph.

    Parameters:
        graph: A DiGraph
        path_list: A list of paths
        path_length: A list of values which are the lengths of the path_list
        weight_avg_list: A list of values which are the weights of the path_list
        delete_entry_node: (Boolean) Delete entry nodes if True
        delete_sink_node : (Boolean) Delete sink nodes if True

    Returns: A DiGraph without paths in path_list exept the best path found.
    """

    # List of index of maximum paths weight (in case of same values)
    max_weight_index = [i for i, j in enumerate(weight_avg_list) if j == max(weight_avg_list)]

    if len(max_weight_index) == 1:
        path_list.pop(max_weight_index[0])
        graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    # If multiple paths has the same max weight, we are going to look for the max length
    else:
        # List of index of maximum paths lengths (in case of same values)
        max_length_index = [i for i, j in enumerate(path_length) if j == max(path_length)]

        if len(max_length_index) == 1:
            path_list.pop(max_length_index[0])
            graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
        # If multiple paths has the same max lengths, we will choose randomly
        else:
            path_list.pop(randint(0, len(path_list) - 1))
            graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph

def path_average_weight(graph, path):
    """
    Compute the mean weight of a path.

    Parameters:
        graph: A DiGraph
        path: A list of paths

    Returns: (Float) The mean weight of the path "path".
    """

    weight = 0

    for i in range(len(path) - 1):
        weight += graph.get_edge_data(path[i], path[i + 1])["weight"]

    return weight / (len(path) - 1)

def solve_bubble(graph, ancestor_node, descendant_node):
    """
    Remove bubbles between ancestor_node and descendant_node inside the graph.

    Parameters:
        graph: A DiGraph
        ancestor_node: A node from the graph
        descendant_node: A node from the graph

    Returns: A DiGraph without bubbles between ancestor_node and descendant_node.
    """

    # Get every possible contig from ancestor_node and descendant_node
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))

    if len(path_list) != 0:
        length_list = []
        weight_list = []

        # Compute the length and mean weight of each contig
        for path in path_list:
            length_list.append(len(path))
            weight_list.append(path_average_weight(graph, path))

        graph = select_best_path(graph, path_list, length_list, weight_list)

    return graph

def simplify_bubbles(graph):
    """
    Find and remove bubbles from the graph.

    Parameters:
        graph: A DiGraph

    Returns: A DiGraph without bubbles.
    """

    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            graph = solve_bubble(graph, starting_node, ending_node)

    return graph


def solve_entry_tips(graph, starting_nodes):
    """
    Find and remove entry tips from the graph.

    Parameters:
        graph: A DiGraph
        starting_nodes: A list of starting nodes from the graph

    Returns: A DiGraph without entry tips.
    """

    descendants_nodes = []

    # Get all descendant nodes from a node in starting_nodes
    for node in starting_nodes:
        descendants_nodes.append(list(descendants(graph, node)))

    # Get nearest shared descendants
    shared_descendants = list(set.intersection(*map(set, descendants_nodes)))

    paths = []
    weights = []

    if len(shared_descendants) != 0:
        for node in starting_nodes:
            # Get the shortest path from the node to the first nearest shared descendant
            path = nx.shortest_path(graph, node, shared_descendants[0])
            paths.append(path)

            weight = 0
            for i in range(len(path) - 1):
                # Compute the weight for each path found
                weight += graph.get_edge_data(path[i], path[i + 1])["weight"]

            weights.append(weight)

        # Remove the path which have the max weight
        # (The first occurence wille be removed if multiple paths has the same max weight)
        paths.pop(weights.index(max(weights)))

        # Remove the others paths from the graph
        graph = remove_paths(graph, paths, True, False)

    return graph

def solve_out_tips(graph, ending_nodes):
    """
    Find and remove out tips from the graph.

    Parameters:
        graph: A DiGraph
        ending_nodes: A list of sink nodes from the graph

    Returns: A DiGraph without out tips.
    """

    predecessors = []

    # Get all ancestors nodes from a node in ending_nodes
    for node in ending_nodes:
        predecessors.append(list(ancestors(graph, node)))

    # Get nearest shared ancestors
    shared_predecessors = list(set.intersection(*map(set, predecessors)))

    paths = []
    weights = []

    if len(shared_predecessors) != 0:
        for node in ending_nodes:
            # Get the shortest path from the node to the last nearest shared descendant
            path = nx.shortest_path(graph, shared_predecessors[-1], node)
            paths.append(path)

            weight = 0
            for i in range(len(path) - 1):
                # Compute the weight for each path found
                weight += graph.get_edge_data(path[i], path[i + 1])["weight"]

            weights.append(weight)

    # Remove the path which have the max weight
    # (The first occurence wille be removed if multiple paths has the same max weight)
    paths.pop(weights.index(max(weights)))

    # Remove the others paths from the graph
    graph = remove_paths(graph, paths, False, True)

    return graph

def get_starting_nodes(graph):
    """
    Get the starting nodes from a graph (where the in_degree equals 0).

    Parameters:
        graph: A DiGraph

    Returns: A list of starting nodes.
    """

    return [n for n, d in graph.in_degree() if d == 0]

def get_sink_nodes(graph):
    """
    Get the sink nodes from a graph (where the out_degree equals 0).

    Parameters:
        graph: A DiGraph

    Returns: A list of sink nodes.
    """

    return [n for n, d in graph.out_degree() if d == 0]

def get_contigs(graph, starting_nodes, ending_nodes):
    """
    Get all possible paths between every nodes in starting_nodes and ending_nodes.

    Parameters:
        graph: A DiGraph
        starting_nodes: A list of starting nodes from the graph
        ending_nodes: A list of sink nodes from the graph

    Returns: A list of tuples with every possible paths with their relative length.
    """

    list = []

    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            # Get all possible paths
            paths = nx.all_simple_paths(graph, starting_node, ending_node)

            for _, path in enumerate(paths):
                if isinstance(starting_node, int):
                    list.append((path, len(path)))
                # If nodes are instance of String
                # then concatenate every path character to make a string path
                else:
                    real_path = path[0]
                    for i in range(1, len(path) - 2):
                        real_path += path[i][1]
                    real_path += path[len(path) - 1]

                    list.append((real_path, len(real_path)))

    return list


def save_contigs(contigs_list, output_file):
    """
    Save contigs in a file.

    Parameters:
        contigs_list: A list of contigs
        output_file: (String) Path to the output file
    """

    # Open the file in "write" mode
    file = open(output_file, "w")

    for index, tuple in enumerate(contigs_list):
        file.write(">contig_" + str(index) + " len=" + str(tuple[1]) + "\n")
        file.write(fill(tuple[0]) + "\n")

    file.close()

def fill(text, width = 80):
    """
    Split text with a line return to respect fasta format
    """

    return os.linesep.join(text[i : i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """
    Draw the graph
    """

    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]

    # Draw the graph with networkx
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')

    # save image
    plt.savefig(graphimg_file)

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    fastq_file = args.fastq_file
    kmer_size = args.kmer_size
    output_file = args.output_file

    dict = build_kmer_dict(fastq_file, kmer_size)
    graph = build_graph(dict)

    graph_without_bubbles = simplify_bubbles(graph)

    starting_nodes = get_starting_nodes(graph_without_bubbles)
    graph_without_entry_tips = solve_entry_tips(graph_without_bubbles, starting_nodes)

    ending_nodes = get_sink_nodes(graph_without_entry_tips)
    graph_without_out_tips = solve_entry_tips(graph_without_entry_tips, ending_nodes)

    draw_graph(graph_without_out_tips, "Graph.jpg")

    starting_nodes = get_starting_nodes(graph_without_out_tips)
    ending_nodes = get_sink_nodes(graph_without_out_tips)
    save_contigs(get_contigs(graph_without_out_tips, starting_nodes, ending_nodes), output_file)

if __name__ == '__main__':
    main()
