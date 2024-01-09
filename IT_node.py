# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 13:54:56 2021

@author: Ali Jazayeri
"""


from IT_nil import nil

class IT_node(nil):
    def __init__(self, interval,node_id):
        nil. __init__(self)
        self.interval = interval
        self.key = interval[0]
        self.left = nil()
        self.right = nil()
        self.maximum = interval[1]
        self.id = node_id
        self.color = 'red'

    def find_max(self):
        max_list = [self.interval[1]]
        if type(self.left) != nil:
            max_list.append(self.left.maximum)
        if type(self.right) != nil:
            max_list.append(self.right.maximum)
        return max(max_list)