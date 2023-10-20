#!/usr/bin/python3
#
from functools import reduce
from uac.edge import Edge


class Element:
    num = 0

    def __init__(self, t, n, reg=0):
        self.type = t
        self.index = Element.num
        self.n = list(map(int, n))
        self.nset = frozenset(self.n)
        self.isolated = False
        self.region(reg)
        self.ctr = None
        self.orig_n = list(self.n)
        self.them = list()
        Element.num += 1

    def __eq__(self, other):
        return isinstance(other, Element) and self.index == other.index

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.index)

    def __repr__(self):
        return "Element(" + self.type + " " + " ".join(map(str, self.n)) + " " + str(self.reg) + ")"

    def region(self, r):
        self.reg = int(float(r))

    def edges(self):
        return [Edge(self.n[i], self.n[(i + 1) % len(self.n)]) for i in range(len(self.n))]

    def fibdir(self, f):
        self.fibdir = f

    def centre(self, pts):
        if self.ctr is None:
            self.ctr = pts[self.n].sum(0) / len(self.n)
        return self.ctr

    def node2ctr(self, data):
        """
        convert node data to center
        """
        avg = reduce(lambda x, y: x + data[y], self.n, 0)
        return avg / len(self.n)

    def has_edge(self, e):
        """
        does element have this edge
        e : edge
        """
        return e.n[0] in self.nset and e.n[1] in self.nset

    def touch_edge(self, e):
        """
        one node in edge
        """
        return (e.n[0] in self.nset) != (e.n[1] in self.nset)

    def neighbour(self, other):
        """
        does other share an edge with self?
        """
        return len(self.nset & other.nset) >= 2

    def is_anti(self, other):
        return other in self.them
