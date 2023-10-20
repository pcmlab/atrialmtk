#!/usr/bin/python3
#


class Line:
    """
    lists of connected edges
    """

    def __init__(self, e=None):
        self.eles = set()  # elements with a node in line
        self.edge = set()  # edges forming line
        self.nodes = set()  # nodes forming line
        if e:
            self.edge.add(e)
            self.nodes = set(e.n)

    def add_edge(self, e):
        self.edge.add(e)
        self.nodes |= set(e.n)

    def join(self, other):
        self.edge |= other.edge
        self.nodes |= other.nodes
