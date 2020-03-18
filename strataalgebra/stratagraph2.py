"""This file is for a a future update...."""
from __future__ import absolute_import

from sage.all import Graph, InfinitePolynomialRing, PolynomialRing
from sage.rings.integer_ring import ZZ


class StrataGraph2(object):
    KappaRing = InfinitePolynomialRing(ZZ, "K")
    K = KappaRing.gen()

    PsiRing = PolynomialRing(ZZ, "psi")
    psi = PsiRing.gen()

    def __init__(self, strataG):
        #make the graph

        G = Graph()

        for v in range(1,strataG.num_vertices()):
            G.add_vertex(v)
            for expon, coef in strataG.M[v,0].dict().items():
                if expon[0] == 1:
                    genus = coef
                else:
                    pass

        self.decorations = dict()

        dec_items = list(self.decorations.items())
        dec_items.sort(lambda x: x[1])

        parts = []
        prev = None
        dec_list = []

        for a in dec_items:
            if a != prev:
                dec_list.append(a[1])
                parts.append(new_part)
                new_part = [a[0]]
            else:
                new_part.append(a[0])

        self.dec_list = tuple(dec_list)

        self.graph, cert = graphUncan.canonical_labeling(parts).copy(immutable = True)

        self.parts = tuple(( tuple([cert[i] for i in part_j].sort) for part_j in parts))

    def __eq__(self,other):
        return self.parts == other.parts and self.dec_list == other.dec_list and self.graph == other.graph

    def __hash__(self):
        return hash(self.parts, self.graph, self.dec_list)
