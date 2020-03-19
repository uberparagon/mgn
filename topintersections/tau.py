"""
Computes the Witten tau function, i.e. the top intersections of psi classes.

Uses the genus recursion of Liu and Xu from [LX07].
"""
from __future__ import absolute_import

from sage.all import Integer, Rational, prod, binomial, multinomial, cartesian_product
from numpy import zeros, array, hstack
from .remember import remember_convert_args


def ee(k, n):
    """
    Returns the k-th n-dimensional "basis vector". 
    """
    v = zeros(n, dtype = Integer)
    v[k] = 1
    return v


#@breadth_first_tree
def tau_no_g(a):
    """
    Like tau, but computes the genus that makes it non-zero, if possible.
    """
    h = 3 + sum([(i-1)*ai for i,ai in enumerate(a)])
    if h % 3 != 0:
        return 0
    else:
        return tau(h/3, a)
        
#@breadth_first_tree
#@remember
#stored_values = dict()

master_table = dict() 

@remember_convert_args(lambda g,a: (g, tuple(a)), master_table) #shares the master_table with intersect_monomial_with_data
def tau(g,a):
    """
    INPUT:
     - g -- The genus.
     - a -- A vector.  The first entry is the number of psi classes with exponent 0, the second in the number with exponent 1, etc.
     
     Uses the recursion of Liu and Xu [LX07].
    """
    
    n = sum(a)
    deg = sum([ai*i for i, ai in enumerate(a)])
    if deg != 3*g-3+n:
        return 0
    if g == 0:
        return multinomial(get_exp_list(a))
    if g == 1 and deg == 1:
        return Rational((1,24))
    if n > 0:
        if a[1] !=0: #dialoton
            return (2*g - 3 + n) * tau(g,a - ee(1,len(a)))
        if a[0] !=0: #string
            return sum( [(ai) * tau(g,a - ee(0, len(a)) + ee(i, len(a)) - ee(i+1,len(a))) for i, ai in enumerate(a[1:]) if ai != 0 ] )


    d = Integer(a.nonzero()[0][0])
    last_index = a.nonzero()[0][-1]
    a1 = array(hstack( (a[0:(last_index+1)],array([0])))) 
    k = len(a1)
    
    a1 -= ee(d, k)    

    ans = Rational((2*d + 3,12)) * tau(g-1, 4*ee(0, k) + ee(d + 1, k) + a1) \
        - Rational((2*g + n -1,6)) * tau(g-1, 3*ee(0,len(a)) + a) \
        + sum(( (binomial_product(a1,b)) * \
            (
                  (2*d + 3) *
                  tau_no_g(2*ee(0,k) + ee(d+1, k) + b) *
                  tau_no_g(2*ee(0,k) + a1 - b) 
                - (2*g + n - 1) *
                  tau_no_g(ee(0,k) + ee(d,k) + b) *
                  tau_no_g(2*ee(0,k) + a1 - b)
            ) for b in splittings(a1)))
            
    return Rational((1,(2*g+n-1)*(2*g+n-2)))*ans

    
def binomial_product(v1,v2):
    """
    See implementation.
    """
    #return prod([binomial(a,b) for a,b in zip(v1, v2)])
    #had to change after update broke it
    return prod([binomial(int(a),int(b)) for a,b in zip(v1, v2)])
    
def splittings(a):
    """
    Important for tau.
    """
    #return [array(c) for c in CartesianProduct(*[range(i+1) for i in a])]
    #sage update broke this
    return [array(c) for c in cartesian_product([list(range(i+1)) for i in a])]
    
def get_exp_list(a):
    """
    Turns a Witten tau list into my traditional list of the exponents of the psis.
    """
    l = []
    for i, ai in enumerate(a):
        l +=[i]*ai
    return l
    
def psi_intersect(g, n, a):
    """
    Here ``a`` is a traditional (to me) list of the exponents of the psis.  This will convert it into a Witten tau list and call the tau function, and return the answer.
    """
    if 3*g - 3 + n != sum(a):
        return 0
     
    d = from_exp_list_to_tau(g,n,a)
     
    return tau(g, array(d, dtype = Integer))
    
def from_exp_list_to_tau(g,n,a):
    """
    Converts a list of psi exponents into a Witten tau list.
    """
    d = [0] * (max(a) +1)
    for ai in a:
        d[ai] += 1
    d[0] = n - len([ai for ai in a if ai != 0])
    return d


    
    
    
