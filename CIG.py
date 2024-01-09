# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 14:03:54 2021

@author: Ali Jazayeri
"""


from collections import defaultdict

import CIG_edge as cig_edge
import CIG_node as cig_node

class cig_network():
    def __init__(self, net_id, directed = False):


        self.net_id = net_id
        self.directed = directed
        # self.edge_counter = -1
        self.nodes = {}
        self.node_lbls = defaultdict(set)
        self.edge_lbls = defaultdict(set)


    # def get_num_nodes(self):
    #     return len(self.nodes)

    def add_node(self, node, lbl):
        self.nodes[node] = cig_node.Node(node, lbl)
        self.node_lbls[lbl].add(node)

    def add_edge(self, frm, to, lbl):
        # self.edge_counter += 1

        self.nodes[frm].add_edge(frm, to, lbl)
        self.edge_lbls[lbl].add((frm, to))

    # def is_neighbor(self, node_id_1, node_id_2):
    #     flag = False

    #     for n in self.nodes[node_id_1].edges.values():
    #         if n.node_id_to == node_id_2:
    #             return True
    #     return flag

    def print_cig(self):
        for nodeix, node in self.nodes.items():
            for edgeix, edge in node.edges.items():
                print(str(edge.frm)+", "+ str(edge.to)+", "+ str(node.lbl)+", "+ str(edge.lbl)+", "+ str(self.nodes[edge.to].lbl))

    def generate_cl(self):
        """
        In this function we create a canonical label for the graph.
        In the future, I need to improve this function by using some of the well-known libraries such as pynauty
        I just create something simple, which is definitely not error-free

        Returns
        -------
        cl : str
            It is an string representing the canonical label of the network.

        """
        cl_components = []
        for nodeix, node in self.nodes.items():
            temp_comp = str(node.lbl) + str(node.len_of_neighbors())
            for edgeix, edge in node.edges.items():
                cl_components.append(temp_comp + "-" +str(edge.lbl))
        return "_".join(sorted(cl_components))
