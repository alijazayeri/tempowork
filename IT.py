# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 13:56:00 2021

@author: Ali Jazayeri
"""


from IT_nil import nil
from IT_node import IT_node

class Interval_Tree(object):
    def __init__(self,root):
        self.root = root
        self.root.color = "black"
        self.root.parent = nil()

    def update_max(self, node):
        if node.maximum != node.find_max():
            node.maximum = node.find_max()
        if type(node.parent) != nil:
            self.update_max(node.parent)

    def Left_Rotate(self, x):
        y = x.right
        x.right = y.left
        if type(y.left) != nil:
            y.left.parent = x
        y.parent = x.parent
        if type(x.parent) == nil:
            self.root = y
        else:
            if x == x.parent.left:
                x.parent.left = y
            else:
                x.parent.right = y
        y.left = x
        x.parent = y

        x.maximum = x.find_max()
        if type(y.left) != nil:
            y.left.maximum = y.left.find_max()
        if type(y.right) != nil:
            y.right.maximum = y.right.find_max()
        y.maximum = y.find_max()
        if type(x.parent) != nil:
            self.update_max(x)
        if type(y.parent) != nil:
            self.update_max(y)


    def Right_Rotate(self, x):
        y = x.left
        x.left = y.right
        if type(y.right) != nil:
            y.right.parent = x
        y.parent = x.parent
        if type(x.parent) == nil:
            self.root = y
        else:
            if x == x.parent.right:
                x.parent.right = y
            else:
                x.parent.left = y
        y.right = x
        x.parent = y

        x.maximum = x.find_max()
        if type(y.left) != nil:
            y.left.maximum = y.left.find_max()
        if type(y.right) != nil:
            y.right.maximum = y.right.find_max()
        y.maximum = y.find_max()
        if type(x.parent) != nil:
            self.update_max(x)
        if type(y.parent) != nil:
            self.update_max(y)

    def Fix_Insertion(self, node):
        while node.parent.color == "red" and node != self.root:

            if node.parent == node.parent.parent.left:
                y = node.parent.parent.right
                if y.color == "red":
                    node.parent.color = "black"
                    node.parent.parent.right.color = "black"
                    node.parent.parent.color = "red"
                    node = node.parent.parent
                else:
                    if node == node.parent.right:
                        node = node.parent
                        self.Left_Rotate(node)
                    node.parent.color = "black"
                    node.parent.parent.color = "red"
                    self.Right_Rotate(node.parent.parent)
            else:
                y = node.parent.parent.left
                if y.color == "red":
                    node.parent.color = "black"
                    node.parent.parent.left.color = "black"
                    node.parent.parent.color = "red"
                    node = node.parent.parent
                else:
                    if node == node.parent.left:
                        node = node.parent
                        self.Right_Rotate(node)
                    node.parent.color = "black"
                    node.parent.parent.color = "red"
                    self.Left_Rotate(node.parent.parent)

            if type(node.left) != nil:
                node.left.maximum = node.left.find_max()
            if type(node.right) != nil:
                node.right.maximum = node.right.find_max()
            node.maximum = node.find_max()
            if type(node.parent) != nil:
                self.update_max(node)

        if type(node.left) != nil:
            node.left.maximum = node.left.find_max()
        if type(node.right) != nil:
            node.right.maximum = node.right.find_max()
        node.maximum = node.find_max()
        if type(node.parent) != nil:
            self.update_max(node)

        self.root.color = "black"

    def Add_Node(self, node):

        y = nil()
        x = self.root
        while type(x) != nil:
            y = x
            if node.key < x.key:
                x = x.left
            else:
                x = x.right
        node.parent = y
        if type(y) == nil:
            self.root = node
        elif node.key < node.parent.key:
            node.parent.left = node
        else:
            node.parent.right = node
        node.left = nil()
        node.right = nil()

        node.maximum = node.find_max()
        if type(node.parent) != nil:
            if type(node.parent.left) != nil:
                node.parent.left.maximum = node.parent.left.find_max()
            if type(node.parent.right) != nil:
                node.parent.right.maximum = node.parent.right.find_max()
            self.update_max(node)
        if node.parent == None:
            node.color = 'red'
            return
        if node.parent.parent == None:
            return
        self.Fix_Insertion(node)

    def do_Overlap(self,i,j):
        if i[0] <= j[1] and j[0] <= i[1]:
            return True
        else:
            return False

    def interval_overlap_all(self,root, i):
        lst = []
        if self.do_Overlap(i,root.interval):
            lst.append((root.id, root.interval))
        if type(root.left) != nil and root.left.maximum >= i[0]:
            lst.extend(self.interval_overlap_all(root.left, i))
        if type(root.right) != nil and root.interval[0] <= i[1] and root.right.maximum >= i[0]:
            lst.extend(self.interval_overlap_all(root.right, i))
        return lst