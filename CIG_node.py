# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 14:02:49 2021

@author: Ali Jazayeri
"""


import CIG_edge as cig_edge

class Node():
    def __init__(self, node, lbl):
        self.node = node
        self.lbl = lbl
        self.edges = {}

    def add_edge(self, frm, to, lbl):
        self.edges[to] = cig_edge.Edge(frm, to, lbl)

    def len_of_neighbors(self):
        return len(self.edges)