#!/usr/bin/python3
#
import numpy as np


class Edge:
    def __init__(self, n0, n1):
        self.n = (n0, n1) if n0 < n1 else (n1, n0)
        self.ele = None

    def __eq__(self, other):
        return isinstance(other, Edge) and self.n == other.n

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.n)

    def edge(self, pts):
        e = pts[self.n[1]] - pts[self.n[0]]
        return e / np.linalg.norm(e)

    def elem(self, a):
        self.ele = a

    def antiele(self, a):
        self.anti = a

    def other(self, n):
        return self.n[n == self.n[0]]

    def attached(self, other):
        """
        are edges connected?
        """
        return self.n[0] in other.n or self.n[1] in other.n

    def vtx(self, n):
        """
        is node n in edge?
        """
        return n == self.n[0] or n == self.n[1]
