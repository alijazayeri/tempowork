# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 14:14:46 2021

@author: Ali Jazayeri
"""


import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd



from src.CIG import cig_network as cig
from src.IT_nil import nil
from src.IT_node import IT_node
from src.IT import Interval_Tree




class frequent():
    def __init__(self, file_name, support, isomorphism = 'e',
                 disc_threshold = 0.05, binning=False, nbins = 10):

        # self.__net_cntr = -1 # network counter to count and id networks in the database
        self.__networks = defaultdict(cig) # a dictionary of cig networks to represent the database
        self.__frequent_nodes = defaultdict(set) # a dictioanry to keep the frequent nodes
        self.__frequent_edges = defaultdict(set) # a dictioanry to keep the frequent edges
        self.__tw_cig_id_map = dict() # we keep track of new ids considered for nodes in CIGs relative to temporal networks
        self.__tw_node_lbl_map = dict() # this dict maps the node labels of the original temporal networks to cig labels
        self.__tw_edge_lbl_map = dict() # this dict maps the edge labels of the original temporal networks to cig labels

        self.__subnk_actual_node_id_map = defaultdict(lambda:defaultdict(list))
        self.__subnk_actual_frequent_node_id_map = defaultdict(lambda:defaultdict(list))

        self.__cl_dict = dict() # this is a list containing the canonical labeling of frequent patterns


        self.__find_evolving = False



        self.__binning = binning
        self.__nbins = nbins


        self.min_supp = support # it is the min frequency that the frequency of subgraphs are compared to.

        self.directed = False

        self.frequent_cntr = 0
        self.__disc_thresh = disc_threshold
        self.__precision = 0.0001
        self.__infinity = np.inf

        self.frequent_patterns = []

        self.file_name = file_name

        if isomorphism == 'e':
            # ----> EXACT Iso
            self.__create_cig_ds("e")
        elif isomorphism == 'i':
            # # # ----> INEXACT Iso
            self.__create_cig_ds("i")
        elif isomorphism == 'es':
            # # ----> SEQUENCE-PRESEREVED Iso EXACT
            self.__create_cig_ds_sequence_preserved("e")
        elif isomorphism == 'is':
            # # ----> SEQUENCE-PRESEREVED Iso INEXACT
            self.__create_cig_ds_sequence_preserved("i")


        self.__reorder_labeling()
        self.__find_frequent_edges()
        self.__run_tw()
        # self.supgraphs_support = self.__subnk_actual_node_id_map


    def __get_edge_elements(self, edge):
        """
        this function gets split line of data

        Parameters
        ----------
        edge : list
            this list is a line of data representing an edge.

        Returns
        -------
        v1 : str
            node 1 id.
        v2 : str
            node 2 id.
        a1 : int
            node 1 label.
        el : int
            edge label.
        a2 : int
            node 2 label.
        st : int
            starting point of the continuous edge.
        dt : int
            duration of the continuous edge.

        """

        v1 = edge[1]
        v2 = edge[2]
        a1 = int(edge[3])
        el = int(edge[4])
        a2 = int(edge[5])
        st = int(edge[6])
        dt = int(edge[7])
        return (v1, v2, a1, el, a2, st , dt)

    def __collect_elements_data(self, elements = ["dt"]):
        """
        This function receives a list of edge elements that user wants to collect to discretize.
        In other words, this functio just iterates over the dataset and
        collects the raw data related to the edge elements that user wants to discretize.

        Parameters
        ----------
        elements : list, default: durations
            The default is ['el'].
            It can be any combination of edge elements, limited to :
                "a1", "el", "a2", "st" , "dt"

        Returns
        -------
        elements_dict: a dictionary, keys: elements, values: lists of values in dataset for elements

        """
        edge_elements = ["v1", "v2", "a1", "el", "a2", "st" , "dt"]
        elements_dict = defaultdict(list)
        with open(self.file_name, "r") as f:
            for line in f:
                if line[0] != 't':
                    line_edge_elements = self.__get_edge_elements(line.strip().split(' '))
                    for e in elements:
                        elements_dict[e].append(line_edge_elements[edge_elements.index(e)])
        return elements_dict

    def __discretize(self, elements_dict):
        """
        This function discretize the raw data values for elements in the dataset.
        We use a pretty simple discretization method.

        ? I might want to elaborate this function further later.

        Parameters
        ----------
        elements_dict: dictionary, keys: elements, values: lists of values in dataset for elements

        Returns
        -------
        disctrization_trees: dictionary, keys: elements, values: disctrization interval trees for elements

        """
        disctrization_trees = dict()
        for element in elements_dict:
            elements_dict[element] = np.array(sorted(elements_dict[element]))
            # we find the points on which the discretization should be performed.
            if self.__binning == 'frequency':
                equal_freq = np.interp(np.linspace(0, len(elements_dict[element]), self.__nbins + 1), np.arange(len(elements_dict[element])),np.sort(elements_dict[element]))
                n, disctrization_edges, patches = plt.hist(elements_dict[element], equal_freq, edgecolor='black')
            elif self.__binning == 'width':
                disctrization_edges = list(list(pd.cut(elements_dict[element], bins=self.__nbins, labels=False, retbins=True))[1])
                disctrization_edges[0] = min(elements_dict[element]) - self.__precision
            else:
                disctrization_edges = elements_dict[element][1:][np.diff(elements_dict[element]) > self.__disc_thresh * elements_dict[element][:-1]]
            # length of discret intervals
            disctrization_cnt = len(disctrization_edges)
            for ix, disctrization_edge in enumerate(disctrization_edges):
                if ix == 0:
                    disctrization_trees[element] = Interval_Tree(IT_node([0, disctrization_edge - self.__precision/100], ix-1))
                    if disctrization_cnt > 1:
                        disctrization_trees[element].Add_Node(IT_node([disctrization_edge, disctrization_edges[ix+1] - self.__precision/100], ix))
                elif ix == disctrization_cnt - 1:
                    disctrization_trees[element].Add_Node(IT_node([disctrization_edge, max(elements_dict[element])*(1+self.__disc_thresh)], ix))
                else:
                    disctrization_trees[element].Add_Node(IT_node([disctrization_edge, disctrization_edges[ix+1] - self.__precision/100], ix))
            if disctrization_cnt == 1:
                disctrization_trees[element].Add_Node(IT_node([disctrization_edges[0], max(disctrization_edges)*(1+self.__disc_thresh)], ix+1))
        return disctrization_trees

    def __create_cig_ds(self, iso_type, elements = None):
        """
        This function receives the isomorphism type from user, and create a dataset of cig networks

        Parameters
        ----------
        iso_type : string
            'e': exact isomorphism
            'i': inexact isomorphism

        Returns
        -------
        None.
        This function updates the self.__networks data structure

        """
        if iso_type == 'i':
            if elements:
                elements_dict = self.__collect_elements_data(elements)
            else:
                elements_dict = self.__collect_elements_data() # if no elements provided by user, just performs discretization over the durations
            disctrization_trees = self.__discretize(elements_dict)
        with open(self.file_name, "r") as f:
            net_cntr = -1
            for line in f:
                if line[0] == "t":
                    v_cnt = -1 # initialize vertex identifier for each network
                    net_cntr += 1
                    self.__networks[net_cntr] = cig(net_cntr, self.directed)
                    IT_dict = dict() # a dictionary for insertion and search of interval trees for the current cig network
                    self.__tw_cig_id_map[net_cntr] = dict()
                elif line[0] == "e":
                    v1, v2, a1, el, a2, st , dt = self.__get_edge_elements(line.strip().split(" "))
                    edge_interval = [st, st + dt]
                    if iso_type == "i":
                        temp_dt = disctrization_trees["dt"].interval_overlap_all(disctrization_trees["dt"].root, [dt, dt + self.__precision])
                        dt = int(np.mean(temp_dt[0][1]))

                    # adding new vertex to the cig
                    v_cnt += 1
                    v_id = v_cnt
                    if not self.directed:
                        a1_temp = min(a1,a2)
                        a2_temp = max(a1,a2)
                        self.__networks[net_cntr].add_node(v_cnt, str(a1_temp) +"-"+ str(el) +"-"+ str(dt) +"-"+ str(a2_temp))
                    else:
                        self.__networks[net_cntr].add_node(v_cnt, str(a1) +"-"+ str(el) +"-"+ str(dt) +"-"+ str(a2))
                    self.__tw_cig_id_map[net_cntr][v_id] = (v1, v2)
                    if v1 in IT_dict: # adding potential dummy edges to the v1_id in cig
                        for v in IT_dict[v1].interval_overlap_all(IT_dict[v1].root, edge_interval):
                            delta = st - v[1][0]
                            if self.directed:
                                if v1 == self.__tw_cig_id_map[net_cntr][v[0]][0]:
                                    direction = '11'
                                else:
                                    direction = '21'
                            else:
                                direction = '00'
                            self.__networks[net_cntr].add_edge(v[0], v_id, str(delta)+"-"+direction)
                        IT_dict[v1].Add_Node(IT_node(edge_interval, v_id))
                    else: # adding interval tree for v1 if it does not have any
                        IT_dict[v1] = Interval_Tree(IT_node(edge_interval, v_id))
                    if v2 in IT_dict: # adding potential dummy edges to the v2_id in cig
                        for v in IT_dict[v2].interval_overlap_all(IT_dict[v2].root, edge_interval):
                            delta = st - v[1][0]
                            if self.directed:
                                if v2 == self.__tw_cig_id_map[net_cntr][v[0]][0]:
                                    direction = '12'
                                else:
                                    direction = '22'
                            else:
                                direction = '00'
                            self.__networks[net_cntr].add_edge(v[0], v_id, str(delta)+"-"+direction)
                        IT_dict[v2].Add_Node(IT_node(edge_interval, v_id))
                    else: # adding interval tree for v2 if it does not have any
                        IT_dict[v2] = Interval_Tree(IT_node(edge_interval, v_id))
                else:
                    print("input data is not correctly formatted!")

    def __create_cig_ds_sequence_preserved(self, iso_type, elements = None):
        """
        This function receives the isomorphism type from user, and create a dataset of cig networks

        Parameters
        ----------
        iso_type : string
            'e': exact isomorphism
            'i': inexact isomorphism

        Returns
        -------
        None.
        This function updates the self.__networks data structure

        """
        if iso_type == 'i':
            if elements:
                elements_dict = self.__collect_elements_data(elements)
            else:
                elements_dict = self.__collect_elements_data() # if no elements provided by user, just performs discretization over the durations
            disctrization_trees = self.__discretize(elements_dict)
        with open(self.file_name, "r") as f:
            net_cntr = -1
            for line in f:
                if line[0] == "t":
                    v_cnt = -1 # initialize vertex identifier for each network
                    net_cntr += 1
                    self.__networks[net_cntr] = cig(net_cntr, self.directed)
                    IT_dict = dict() # a dictionary for insertion and search of interval trees for the current cig network
                    self.__tw_cig_id_map[net_cntr] = dict()
                elif line[0] == "e":
                    v1, v2, a1, el, a2, st , dt = self.__get_edge_elements(line.strip().split(" "))
                    edge_interval = [st, st + dt]
                    if iso_type == "i":
                        # print([dt, dt + self.__precision])
                        dt = int(np.mean(disctrization_trees["dt"].interval_overlap_all(disctrization_trees["dt"].root, [dt, dt + self.__precision])[0][1]))

                    # adding new vertex to the cig
                    v_cnt += 1
                    v_id = v_cnt
                    if not self.directed:
                        a1_temp = min(a1,a2)
                        a2_temp = max(a1,a2)
                        self.__networks[net_cntr].add_node(v_cnt, str(a1_temp) +"-"+ str(el) +"-"+ str(dt) +"-"+ str(a2_temp))
                    else:
                        self.__networks[net_cntr].add_node(v_cnt, str(a1) +"-"+ str(el) +"-"+ str(dt) +"-"+ str(a2))
                    self.__tw_cig_id_map[net_cntr][v_id] = (v1, v2)

                    delta = 1
                    if v1 in IT_dict: # adding potential dummy edges to the v1_id in cig
                        for v in IT_dict[v1].interval_overlap_all(IT_dict[v1].root, edge_interval):
                            if self.directed:
                                if v1 == self.__tw_cig_id_map[net_cntr][v[0]][0]:
                                    direction = '11'
                                else:
                                    direction = '21'
                            else:
                                direction = '00'
                            self.__networks[net_cntr].add_edge(v[0], v_id, str(delta)+"-"+direction)
                        IT_dict[v1].Add_Node(IT_node(edge_interval, v_id))
                    else: # adding interval tree for v1 if it does not have any
                        IT_dict[v1] = Interval_Tree(IT_node(edge_interval, v_id))
                    if v2 in IT_dict: # adding potential dummy edges to the v2_id in cig
                        for v in IT_dict[v2].interval_overlap_all(IT_dict[v2].root, edge_interval):
                            if self.directed:
                                if v2 == self.__tw_cig_id_map[net_cntr][v[0]][0]:
                                    direction = '12'
                                else:
                                    direction = '22'
                            else:
                                direction = '00'
                            self.__networks[net_cntr].add_edge(v[0], v_id, str(delta)+"-"+direction)
                        IT_dict[v2].Add_Node(IT_node(edge_interval, v_id))
                    else: # adding interval tree for v2 if it does not have any
                        IT_dict[v2] = Interval_Tree(IT_node(edge_interval, v_id))
                else:
                    print("input data is not correctly formatted!")

    def __reorder_labeling(self):
        """
        This function reorders the labels (gives new labels) of network nodes and edges based on their frequencies.

        Returns
        -------
        None.
            It updates the self.__networks data structure

        """
        # two temporary dictionaries of sets to keep the id of networks having a specific lbl
        node_lbl_freq = defaultdict(set)
        edge_lbl_freq = defaultdict(set)


        for net_ix , net in self.__networks.items():
            for _, node in net.nodes.items():
                node_lbl_freq[node.lbl].add(net_ix)
                for _, edge in node.edges.items():
                    edge_lbl_freq[edge.lbl].add(net_ix)



 	# two vectors are created for nodes and edges composed of pairs of node and edge labels and their size as follows:
 	# 		first one is populated with node labels and the size of each label
 	# 		second one is populated with the edge labels and the size of each label
        temp_nodes_vec = list()
        temp_edges_vec = list()
        for lbl, supp in node_lbl_freq.items():
            temp_nodes_vec.append((lbl, len(supp)))
        for lbl, supp in edge_lbl_freq.items():
            temp_edges_vec.append((lbl, len(supp)))

        # sorting list of node and edge labels based on their frequencies
        temp_nodes_vec.sort(key = lambda x:-x[1])
        temp_edges_vec.sort(key = lambda x:-x[1])

 	# two maps are created for mapping previous labels of nodes and edges
 	# with new labels based on the frequency of labels of nodes and edges
        node_map = dict()
        edge_map = dict()

        node_lbl = -1
        for lbl, _ in temp_nodes_vec:
            node_lbl += 1
            node_map[lbl] = node_lbl
            self.__tw_node_lbl_map[node_lbl] = lbl

        edge_lbl = -1
        for lbl, _ in temp_edges_vec:
            edge_lbl += 1
            if type(lbl) == int:
                edge_map[lbl] = edge_lbl + 1000
                self.__tw_edge_lbl_map[edge_lbl + 1000] = lbl
            else:
                edge_map[lbl] = edge_lbl
                self.__tw_edge_lbl_map[edge_lbl] = lbl

 	# temp_networks is created to reorder labels in the self.__networks based on their frequencies
        temp_networks = defaultdict(cig)
        for net_ix, net in self.__networks.items():
            temp_networks[net_ix] = cig(net_ix, self.directed)
            for node_ix, node in net.nodes.items():
                temp_networks[net_ix].add_node(node.node, node_map[node.lbl])
                for edge_ix, edge in node.edges.items():
                    temp_networks[net_ix].add_edge(edge.frm, edge.to, edge_map[edge.lbl])

        self.__networks = temp_networks.copy()


    def __find_frequent_edges(self):
        """
        In the following function, we iterate through all the networks, then nodes and then edges
        and record the number of appearances of each node and edge label in all the networks.
        If the number of appearance is less than the minimum support threshold, we deleted the node and the edge (they are infrequent).
        We know that infrequent nodes and edges cannot be considered later as an edge for frequent subnetwork (Downward Closure Property).
        the main output of this function is two vectors including frequent nodes and edges


        Returns
        -------
        None.
            This function updates two data structures frequent edges and frequent nodes
        """
        for net_ix, net in self.__networks.items():  # iterates over networks composed of relabeled nodes and edges
            for node_ix, node in net.nodes.items():#iterates over nodes of each network
                self.__frequent_nodes[node.lbl].add(net_ix)
                self.__subnk_actual_frequent_node_id_map[node.lbl][net_ix].append(node.node)
                for edge_ix, edge in node.edges.items(): #iterates over all the edges of each node
                    nd_frm_lbl = net.nodes[edge.frm].lbl
                    eg_lbl = edge.lbl
                    nd_to_lbl = net.nodes[edge.to].lbl
                    temp_tw_input_id = dict()
                    temp_edge = tuple()
                    temp_edge = (nd_frm_lbl, eg_lbl, nd_to_lbl)
                    temp_tw_edge_vect = tuple([(-100001,-100000, temp_edge[0], temp_edge[1], temp_edge[2])])
                    temp_tw_input_id[-100001] = edge.frm
                    temp_tw_input_id[-100000] = edge.to
                    temp_tw_input_id[edge.frm] = -1
                    temp_tw_input_id[edge.to] = -1
                    self.__frequent_edges[temp_edge].add(net_ix)
                    self.__subnk_actual_node_id_map[temp_tw_edge_vect][net_ix].append(temp_tw_input_id)
        for edge in list(self.__subnk_actual_node_id_map):
            if len(self.__subnk_actual_node_id_map[edge]) < self.min_supp:
                del self.__subnk_actual_node_id_map[edge]
        # Here, the frequent_edges is finalized with the list of frequent edges which have frequency more than min_support
        for edge in list(self.__frequent_edges):
            if len(self.__frequent_edges[edge]) < self.min_supp:
                del self.__frequent_edges[edge]
        # Here, the frequent_nodes is finalized with the list of frequent nodes which have frequency more than min_support
        for node in list(self.__frequent_nodes):
            if len(self.__frequent_nodes[node]) < self.min_supp:
                del self.__frequent_nodes[node]

    def __run_tw(self):
        """
        #	It is main function of mining frequent subgraphs.
        #	It sends the frequent edges to be mined
        #	Here, having the list of frequent size one subgraphs (edges), we iterate over them
        #	and following a dfs strategy find the frequent larger subgraphs


        Returns
        -------
        None.

        """
        #	I needed to create one frequent_edges_map similar to subnk_actual_node_id_map
        #	These two include the frequent edges, but after each iteration over frequent_edges_map_reduced, we remove the mined
        #	edge from frequent_edges_map, so it is not used in the next mining operations. This removal narrows down the search space

        self.__frequent_edges_map = self.__subnk_actual_node_id_map.copy()
        seed_edges = [i for i in self.__subnk_actual_node_id_map.keys()]
        seed_edges.sort()

        for node in self.__subnk_actual_frequent_node_id_map:
            if len(self.__subnk_actual_frequent_node_id_map[node]) >= self.min_supp:
                one_support = list(self.__subnk_actual_frequent_node_id_map[node].items())[0]
                v1l, el, dt, v2l = self.__tw_node_lbl_map[node].split("-")
                v1id, v2id = self.__tw_cig_id_map[one_support[0]][one_support[1][0]]
                self.frequent_cntr += 1
                first_edge = [(v1id, v2id,v1l, el, v2l, 0, dt)]
                self.frequent_patterns.append(first_edge)


        for edge in seed_edges:
            edge_supp = self.__frequent_edges_map[edge]
            subnk_rmp_map = edge
            one_support = list(self.__subnk_actual_node_id_map[edge].items())[0]
            reconsted, end_time = self.__reconstruct_tw(edge, one_support)
            cl = self.__construct_cig(reconsted).generate_cl()
            if hash(cl) not in self.__cl_dict:

                self.frequent_patterns.append(reconsted)
                self.frequent_cntr += 1
                frequent_cntr = self.frequent_cntr
                self.__cl_dict[hash(cl)] = frequent_cntr
                if self.__find_evolving:
                    self.frequents_nets_range[self.frequent_cntr] = [0,int(self.__tw_edge_lbl_map[edge[0][-2]].split('-')[1])]
                self.__add_forward_edge(edge, edge_supp, subnk_rmp_map, frequent_cntr)
                if not self.__find_evolving:
                    del self.__frequent_edges_map[edge]

    def __is_valid(self, edge_vect, new_edge):
        """
  	   Here, we make sure that n-size child subgraph going to be created from (n-1)-size edge_vect subgraph and
  	   the new_edge is valid following right-most growth strategies provided by gSpan technical report
  	   They are called DFS Code's Neighborhood Restriction rules in the report

        Parameters
        ----------
        edge_vect : tuple of edge(s)
            it is an edge vector.
        new_edge : tuple
            it is an edge.

        Returns
        -------
        Boolean.
            It returns true if the combination of edge vector and new edge is valid, false otherwise.
        """
        is_valid = False
        if edge_vect[-1][0] < edge_vect[-1][1]:
            if new_edge[0] < new_edge[1]:
                if new_edge[0] <= edge_vect[-1][1] and new_edge[1] == edge_vect[-1][1] + 1:
                    is_valid = True
            if new_edge[0] > new_edge[1]:
                if new_edge[0] > edge_vect[-1][0] and new_edge[1] == edge_vect[-1][1]:
                    is_valid = True
        return is_valid

    def __is_minimum_dfs(self, edge_vect, new_edge):
        """
        This function has (checks) all the pruning strategies recommended by gSpan
        to prune duplicate (or non-minimum) subgraphs
        It receives edge_vect (which is an (n-1)-size subgraph) and a new edge (new_edge)
        It returns true if forming the n-size subgraph from edge_vect by adding the new_edge doesn't produce an non-minimum (or duplicate) subgraph
        It does that following three strategies as below. These strategies are coming from gSpan technical report


        Parameters
        ----------
        edge_vect : tuple of edge(s)
            it is an edge vector.
        new_edge : tuple
            it is an edge.

        Returns
        -------
        Boolean.
            It returns true if forming the n-size subgraph from edge_vect by adding the new_edge doesn't produce an non-minimum (or duplicate) subgraph, false otherwise.
        """
        temp_edge_vect = edge_vect + tuple([new_edge])
        is_valid_code = True
        temp_new_edge = tuple()
        if new_edge[2] < new_edge[4]:
            temp_new_edge = (new_edge[2], new_edge[3], new_edge[4])
        else:
            temp_new_edge = (new_edge[4], new_edge[3], new_edge[2])
  	   # 1. the child of edge_vect (here it is called temp_edge_vect) should not include any edge smaller than
  	   # the first edge of the edge_vect
        if new_edge[0] < new_edge[1]:
            if ( (  (temp_edge_vect[0][0] >= new_edge[0]) or
    				(temp_edge_vect[0][0] == new_edge[0] and temp_edge_vect[0][2] <= new_edge[2]) or
    				(temp_edge_vect[0][0] == new_edge[0] and temp_edge_vect[0][2] == new_edge[2] and abs(temp_edge_vect[0][3]) <= abs(new_edge[3])) or
    				(temp_edge_vect[0][0] == new_edge[0] and temp_edge_vect[0][2] == new_edge[2] and abs(temp_edge_vect[0][3]) == abs(new_edge[3])  and temp_edge_vect[0][4] <= new_edge[4]) or
    				(temp_edge_vect[0][1] < new_edge[1])
    				)):
                is_valid_code = True
            else:
                is_valid_code = False
                return is_valid_code
        else:
            if temp_edge_vect[0][1] <= new_edge[0]:
                is_valid_code = True
            else:
                is_valid_code = False
                return is_valid_code


        # If there is an edge connected by the head to node N,
        # It should NOT be smaller than all other edges connected by their head to node N
        for i in range(len(temp_edge_vect)-1):
            if temp_edge_vect[i][1] == new_edge[1]:
                if (
                    (temp_edge_vect[i][2] > new_edge[2]) or
                    (temp_edge_vect[i][2] == new_edge[2] and temp_edge_vect[i][3] > new_edge[3]) or
                    (temp_edge_vect[i][2] == new_edge[2] and temp_edge_vect[i][3] == new_edge[3] and temp_edge_vect[i][4] > new_edge[4])
                        ):
                    is_valid_code = False
                    return is_valid_code

        # If there is an edge connected by tail to (orginated from) node N,
        # It should NOT be smaller than all other edges orginated from node N
        for i in range(len(temp_edge_vect)-1):
            if temp_edge_vect[i][0] == new_edge[0]:
                if (
                    (temp_edge_vect[i][2] > new_edge[2]) or
                    (temp_edge_vect[i][2] == new_edge[2] and temp_edge_vect[i][3] > new_edge[3]) or
                    (temp_edge_vect[i][2] == new_edge[2] and temp_edge_vect[i][3] == new_edge[3] and temp_edge_vect[i][4] > new_edge[4])
                        ):
                    is_valid_code = False
                    return is_valid_code
        return is_valid_code


    def __add_forward_edge(self, edge_vect, edge_vetw_supp, subnk_rmp_map, frequent_cntr):
        """
    	   Here, we try to add forward edges to the nodes of edge_vect located on the right-most path,
    	   starting from right-most node (node n) to the first node (node 0) of right-most path

        Parameters
        ----------
        edge_vect : tuple of edge(s)
            it is an edge vector which is going to be extended by a forward edge.
        edge_vetw_supp : dictionary
            it contains the supporting networks and their address of edge_vect in the dataset.
        temp_subnk_rmp_map : edge vector
            It is the right most path of edge vector.
        frequent_cntr : int
            it is an integer representing the numver of frequent patterns identified so far.

        Returns
        -------
        None.

        """
        temp_subnk_rmp_map = subnk_rmp_map
        for i in range(len(temp_subnk_rmp_map) -1, -2, -1):
            if i > -1:
                l_nd_ix = 0
                r_nd_ix = 1
                r_nd_lbl = 4
                it_node = i
            elif i == -1:
                l_nd_ix = 1
                r_nd_ix = 0
                r_nd_lbl = 2
                it_node = 0
         	  #	We iterate over the frequent edges to check different possibilities of adding a new edge to the
         	  #	edge_vect as forward edges
            for seed, seed_supp in self.__frequent_edges_map.items():

                temp_subnk_actual_node_id_map = defaultdict(lambda:defaultdict(list))
                if subnk_rmp_map[it_node][r_nd_lbl] == seed[0][2]:
                    temp_edge_vect = edge_vect
                    temp_edge = (edge_vect[it_node][r_nd_ix], edge_vect[-1][1] + 1, seed[0][2], seed[0][3], seed[0][4])
                    checks = False

                    if self.__is_valid(temp_edge_vect, temp_edge):
                        if self.__is_minimum_dfs(temp_edge_vect, temp_edge):
                            checks = True
                    if checks:

                        temp_edge_vect += tuple([temp_edge])

                        temp_rmp = tuple() # ? any better way than for loop
                        for forward_ix in range(i+1):
                            temp_rmp += tuple([temp_subnk_rmp_map[forward_ix]])
                        temp_rmp += tuple([temp_edge])
                        for supp_ix, supp in edge_vetw_supp.items():
                            for supp_map in supp:
                                nd_ids = set()
                                for _ , mp in supp_map.items():
                                    nd_ids.add(mp)
                                for freq_edge_supp_map in seed_supp[supp_ix]:
                                    temp_map = supp_map.copy()
                                    if freq_edge_supp_map[seed[0][1]] in nd_ids:
                                        is_in = True
                                    else:
                                        is_in = False
                                    if supp_map[temp_subnk_rmp_map[it_node][r_nd_ix]] == freq_edge_supp_map[seed[0][0]] and  supp_map[temp_subnk_rmp_map[it_node][l_nd_ix]] != freq_edge_supp_map[seed[0][1]] and not is_in:
                                        temp_map[temp_edge[1]] = freq_edge_supp_map[seed[0][1]]
                                        temp_subnk_actual_node_id_map[temp_edge_vect][supp_ix].append(temp_map)

                        if (len(temp_subnk_actual_node_id_map[temp_edge_vect]) >= self.min_supp and temp_edge_vect != edge_vect):

                                one_support = list(temp_subnk_actual_node_id_map[temp_edge_vect].items())[0]
                                reconsted, end_time = self.__reconstruct_tw(temp_edge_vect, one_support)
                                cl = self.__construct_cig(reconsted).generate_cl()
                                if hash(cl) not in self.__cl_dict:
                                    self.__subnk_actual_node_id_map[temp_edge_vect] = temp_subnk_actual_node_id_map[temp_edge_vect]

                                    self.frequent_patterns.append(reconsted)
                                    self.frequent_cntr += 1
                                    frequent_cntr = self.frequent_cntr
                                    self.__cl_dict[hash(cl)] = frequent_cntr
                                    self.__add_forward_edge(temp_edge_vect, temp_subnk_actual_node_id_map[temp_edge_vect], temp_rmp, frequent_cntr)


                elif subnk_rmp_map[it_node][r_nd_lbl] == seed[0][4]:
                    temp_edge_vect = edge_vect
                    temp_edge = (edge_vect[-1][1] + 1, edge_vect[it_node][r_nd_ix], seed[0][2], seed[0][3], seed[0][4])
                    checks = False
                    if self.__is_valid(temp_edge_vect, temp_edge):
                        if self.__is_minimum_dfs(temp_edge_vect, temp_edge):
                            checks = True
                    if checks:
                        temp_edge_vect += tuple([temp_edge])
                        temp_rmp = tuple() # ? any better way than for loop
                        for forward_ix in range(i+1):
                            temp_rmp += tuple([temp_subnk_rmp_map[forward_ix]])
                        temp_rmp += tuple([temp_edge])
                        for supp_ix, supp in edge_vetw_supp.items():
                            for supp_map in supp:
                                nd_ids = set()
                                for map_ix, mp in supp_map.items():
                                    nd_ids.add(mp)
                                for freq_edge_supp_map in seed_supp[supp_ix]:
                                    temp_map = supp_map.copy()
                                    if freq_edge_supp_map[seed[0][0]] in nd_ids:
                                        is_in = True
                                    else:
                                        is_in = False
                                    if supp_map[temp_subnk_rmp_map[it_node][r_nd_ix]] == freq_edge_supp_map[seed[0][1]] and \
                                    supp_map[temp_subnk_rmp_map[it_node][l_nd_ix]] != freq_edge_supp_map[seed[0][0]] and not is_in:
                                        temp_map[temp_subnk_rmp_map[-1][1] + 1] = freq_edge_supp_map[seed[0][0]]
                                        temp_subnk_actual_node_id_map[temp_edge_vect][supp_ix].append(temp_map)
                        # if len(edge_vect) == 1:
                        #     print(edge_vect, seed,"-->",len(temp_subnk_actual_node_id_map[temp_edge_vect]))
                        if (len(temp_subnk_actual_node_id_map[temp_edge_vect]) >= self.min_supp and temp_edge_vect != edge_vect):

                                one_support = list(temp_subnk_actual_node_id_map[temp_edge_vect].items())[0]
                                reconsted, end_time = self.__reconstruct_tw(temp_edge_vect, one_support)
                                cl = self.__construct_cig(reconsted).generate_cl()
                                if hash(cl) not in self.__cl_dict:
                                    self.__subnk_actual_node_id_map[temp_edge_vect] = temp_subnk_actual_node_id_map[temp_edge_vect]
                                    # print("------------")
                                    # print(temp_edge_vect)
                                    # print(reconsted)
                                    self.frequent_patterns.append(reconsted)
                                    self.frequent_cntr += 1
                                    frequent_cntr = self.frequent_cntr
                                    self.__cl_dict[hash(cl)] = frequent_cntr
                                #     self.__add_forward_edge(temp_edge_vect, temp_subnk_actual_node_id_map[temp_edge_vect], temp_rmp, frequent_cntr)
                                    self.__add_forward_edge(temp_edge_vect, temp_subnk_actual_node_id_map[temp_edge_vect], temp_rmp, frequent_cntr)


    def __create_edges_dict(self, edge_vect):
        """


        Parameters
        ----------
        edge_vect : tuple of edge(s)
            it is an edge vector.

        Returns
        -------
        neighbors_dict : dictionary
            a dictionary of nodes:[list of neighbors]
        lbls_dict : dictionary
            a dictionary recording the lbls of nodes and edges

        """
        edges_dict = defaultdict(lambda:defaultdict(list))
        for edge in edge_vect:
            edges_dict[edge[0]]['neighbor'].append(edge[1])
            edges_dict[edge[1]]['neighbor'].append(edge[0])
            edges_dict[edge[0]]['lbl'] = edge[2]
            edges_dict[(edge[0], edge[1])]['lbl'] = edge[3]
            edges_dict[(edge[1], edge[0])]['lbl'] = edge[3]
            edges_dict[(edge[0], edge[1])]['dir'] = 1
            edges_dict[(edge[1], edge[0])]['dir'] = -1
            edges_dict[edge[1]]['lbl'] = edge[4]
        return edges_dict

    def __reconstruct_tw(self, edge_vect, one_support):
        """
        This function receives an edge vector and return the temporal network representing the edge_vector
        The dummy egdes are used to find the delays or differences between pair of original edges connected by each dummy edge.

        Parameters
        ----------
        edge_vect : tuple of edge(s)
            it is an edge vector.
        one_support : dictionary
            it represents one of them supports of edge_vect in one of the networks in ds.

        Returns
        -------
        temporal_network: a list of edges
            it is formatted similar to one of the networks in the dataset

        """
        edges_dict = self.__create_edges_dict(edge_vect)
        sts = {} # starting times of nodes
        temporal_network = [] # temporal network
        source = edge_vect[0][0]
        e1v1id, e1v2id = self.__tw_cig_id_map[one_support[0]][one_support[1][0][source]] # the initial node added to the edge_vect
        sts[source] = 0
        e1v1l,e1l,dt1,e1v2l = self.__tw_node_lbl_map[edges_dict[source]['lbl']].split('-')
        e1v1id, e1v2id = self.__tw_cig_id_map[one_support[0]][one_support[1][0][source]]
        temporal_network.append((e1v1id, e1v2id,e1v1l,e1l,e1v2l,sts[source], dt1))
        stack = [source]  # a list of nodes visited
        end_time = []
        visited = []
        negative_time = 0
        while len(stack) != 0:
            s = stack.pop()
            if s not in visited:
                visited.append(s)
                e1v1id, e1v2id = self.__tw_cig_id_map[one_support[0]][one_support[1][0][s]]
                for neighbor in edges_dict[s]['neighbor']:
                    if neighbor not in visited:
                        e2v1l,e2l,dt2,e2v2l = self.__tw_node_lbl_map[edges_dict[neighbor]['lbl']].split('-')
                        e2v1id, e2v2id = self.__tw_cig_id_map[one_support[0]][one_support[1][0][neighbor]]
                        duration, direction = self.__tw_edge_lbl_map[edges_dict[(s, neighbor)]['lbl']].split('-')
                        sts[neighbor] = sts[s] + edges_dict[(s, neighbor)]['dir'] * int(duration)
                        if sts[neighbor] < 0:
                            negative_time = min(negative_time, sts[neighbor])
                        temp_e2v1id = e2v1id
                        temp_e2v2id = e2v2id
                        temporal_network.append((temp_e2v1id, temp_e2v2id,e2v1l,e2l,e2v2l,sts[neighbor], dt2))
                        end_time.append(sts[neighbor])
                        stack.append(neighbor)
        if negative_time < 0:
            updated_temporal_net = []
            updated_end_time = []
            for ix, tup in enumerate(temporal_network):
                temp_tup = list(tup)
                temp_tup[-2] -= negative_time
                updated_temporal_net.append(temp_tup)
                updated_end_time.append(temp_tup[-2])
            temporal_network = updated_temporal_net
            end_time = updated_end_time
        return sorted(temporal_network, key= lambda x:(x[-2],x[-1],x[2],x[3],x[4])), max(end_time)

    def __construct_cig(self, tw_list):
        """
        This function creates a CIG from a list of temporal edges

        Parameters
        ----------
        tw_list : list
            it is a list of temporal edges representing a temporal network.

        Returns
        -------
        cig_net : CIG network
            It is the CIG representing the temporal network of the tw_list.

        """
        v_id = -1 # initialize vertex identifier for each network
        net = cig(0, self.directed)
        IT_dict = dict() # a dictionary for insertion and search of interval trees for the current cig network
        pairs_dict = dict()
        for line in tw_list:

            v1, v2, a1, el, a2, st , dt = line
            edge_interval = [int(st), int(st) +  int(dt)]
            # adding new vertex to the cig
            v_id += 1
            if not self.directed:
                a1_temp = min(a1,a2)
                a2_temp = max(a1,a2)
                net.add_node(v_id, str(a1_temp) +"-"+ str(el) +"-"+ str(dt) +"-"+ str(a2_temp))
            else:
                net.add_node(v_id, str(a1) +"-"+ str(el) +"-"+ str(dt) +"-"+ str(a2))
            pairs_dict[v_id] = (v1, v2)
            if v1 in IT_dict: # adding potential dummy edges to the v1_id in cig
                for v in IT_dict[v1].interval_overlap_all(IT_dict[v1].root, edge_interval):
                    delta = st - v[1][0]
                    if self.directed:
                        if v1 == pairs_dict[v[0]][0]:
                            direction = '11'
                        else:
                            direction = '21'
                    else:
                        direction = '00'
                    net.add_edge(v[0], v_id, str(delta) + '-' + direction)
                IT_dict[v1].Add_Node(IT_node(edge_interval, v_id))
            else: # adding interval tree for v1 if it does not have any
                IT_dict[v1] = Interval_Tree(IT_node(edge_interval, v_id))
            if v2 in IT_dict: # adding potential dummy edges to the v2_id in cig
                for v in IT_dict[v2].interval_overlap_all(IT_dict[v2].root, edge_interval):
                    delta = st - v[1][0]
                    if self.directed:
                        if v2 == pairs_dict[v[0]][0]:
                            direction = '21'
                        else:
                            direction = '22'
                    else:
                        direction = '00'
                    net.add_edge(v[0], v_id, str(delta) + '-' + direction)
                IT_dict[v2].Add_Node(IT_node(edge_interval, v_id))
            else: # adding interval tree for v2 if it does not have any
                IT_dict[v2] = Interval_Tree(IT_node(edge_interval, v_id))
        return net