"""
The ``tests`` module tests the implementation of the product and the FZ relations by comparing them to the ``topintersections`` code.

You can change which pairs you want to check by modifying the ``g_n_pairs_to_check variable`` in the source code ``tests.py``.

::

    sage: import strataalgebra.tests as tests
    sage: tests.run() # not tested

The comment is there so that the doctests don't try to run it.
"""

from sage.all import QQ, floor, WeightedIntegerVectors, subsets
#from sage.all import *
from strataalgebra import StrataAlgebra

from topintersections import intnum

g_n_pairs_to_check = [
(2,0),
(1,1),
(1,2),
(1,3),
(2,1),
(2,2)
]

def get_all_combis(g,n):
    dim = 3*g-3 + n
    reducible_boundaries = 0
    marks = range(1,n+1)
    if n != 0:
        first_mark_list = [marks.pop()]            
        for g1 in range(0, g + 1):
            for p in subsets(marks):
                r_marks = set(first_mark_list + p)
                if 3*g1 - 3 + len(r_marks) + 1 >= 0 and 3*(g-g1) - 3 + n - len(r_marks) + 1 >= 0:
                    reducible_boundaries+=1  
    
    else: #self.n == 0
        for g1 in range(1, floor(g/2.0)+1):
            reducible_boundaries+=1

    #print "computed red bound"
    indexes = range(1,n+dim+1) + range(n+dim+g+1, n+dim+g+reducible_boundaries + 2)
    codims = [1]*n + range(1,dim+1) + [1]*(reducible_boundaries +1)
    
    for w in WeightedIntegerVectors(dim,codims):
        combi = []
        #print w
        for index, wi in zip(indexes,w):
            combi += [index]*wi
        yield combi

def run():
    """
    run doc
    :return:
    """
    ok = 0
    bad = 0
    for g,n in g_n_pairs_to_check:
        print "Beginning to check g = {0}, n = {1}...".format(g,n)
        s = StrataAlgebra(QQ,g,tuple(range(1,n+1)),False)
        for c in get_all_combis(g,n):
            print c,
            v1 = intnum(g,n,c, confirm = False)
            v2 = s.MgnLb_int(c)
            if v1==v2:
                print " ok."
                ok+=1
            else:
                print
                print " FAILED!"
                bad+=1
                #print c
                print v1, "vs", v2

    print "{0} tests performed, {1} successful, {2} failed".format(ok+bad,ok,bad)
            
class doctestclass(object):
    """
    here is the doc
    """
    pass