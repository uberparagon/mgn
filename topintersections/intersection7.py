"""
This is the file with the main intersection code!
"""
from __future__ import absolute_import, print_function

try:
    from sage.all import *
    from .TautRing3 import *
    from .tau import *
    #from checkin import *
    from .remember import remember_convert_args
    from .push_down import push_down


except ImportError:
    pass

    
class monomial_data(object):
    """
    This class parses a monomial and stores the important data about it.  The data can then be easily passed around.
    """

    def __init__(self, space, monomial=None):
        """
        This constructor initializes the data and parses the monomial.  
        
        If you don't pass in a monomial, then it is assumed that you will call the recieve method manually.  This happens in instances of monomial_data_product.
        """
        self.space = space
        self.psis = []
        self.red_boundaries = []
        self.irr_boundaries = 0
        self.kappas = []
        self.cherns = []
        self.lambdas = []
        if None != monomial:
            self.monomial = monomial
            for tclass, expon in monomial.decompose_monomial():
                self.recieve(tclass, expon)
            self.recieve = None #make sure we don't double populate it  
        else:
            self.monomial = 1
            
    #########################################################################
    #The following methods are used to do the parsing and intialize the data.
    
    def psi(self, tclass, expon):
        self.psis.append((tclass, expon))
        
    def irr(self, tclass, expon):
        self.irr_space = tclass.space
        self.irr_boundaries += expon
        
    def red(self, tclass, expon):
        self.red_boundaries += [tclass]*expon
        
    def kappa(self, tclass, expon):
        self.kappas += [tclass.degree]*expon
        
    def chern(self, tclass, expon):
        self.cherns += [tclass.degree]*expon
        
    def lamb(self, tclass, expon):
        self.lambdas += [tclass.degree]*expon
        
    # This is a class variable (dict) that associates a method to each class.
    route = {   psi_class: psi,
                reducible_boundary: red,
                irreducible_boundary: irr,
                kappa_class: kappa,
                chern_char: chern,
                lambda_class: lamb }
                
    #  end parsing methods
    ##################################################################
    
    def recieve(self, tclass, expon):
        """
        Recieves a TautRingElement (tclass) and an integer exponent, and calls the appropriate method to record the data.
        """
        self.route[tclass.__class__](self, tclass, expon)
        
    def psi_exp(self):
        """
        Returns a list of the exponents of the psis.
        """
        return [t[1] for t in self.psis]
        
    def degree(self):
        """
        Computes the degree of the class represented by this monomial.
        """
        return sum(self.psi_exp()) + sum(self.kappas) + len(self.red_boundaries) + self.irr_boundaries + sum(self.cherns) + sum(self.lambdas)        
    
    def get_chern(self):
        """
        Returns the chern_char class that you should deal with first.
        
        TODO:  decide which one you should do first. for now, it just returns the first one.
        """
        return chern_char(self.space, self.cherns[0])
    
    def get_lambda(self):
        """
        Returns the lambda_class that you should substitute out first, and its exponent.
        
        TODO:  decide which one you should do first.
        """
        return lambda_class(self.space, self.lambdas[0]), self.lambdas.count(self.lambdas[0])
        
    def get_red(self):
        """
        Returns the reducible_boundary class that you should pull back onto first.
        
        Assume it has already been sorted since self.key will call sort(self.irr_boundaries).  See the reducible_boundary class in the file TautRing3.py for the comparison methods.
        """
        return self.red_boundaries[0]
        
    def vanishes_by_degree(self):
        """
        Returns True if the class represented by this monomial does not have degree equal to the space it is on.  Returns False if it does have the right degree.
        """
        return self.degree() != self.space.dimension
        
    def base_case(self):
        """
        Checks if this monomial is covered by some base cases.  If it is, then return the value.  If not, return None.
        
        TODO:  Is this the best way to handle this, or should I have a dict or something?
        """
        if self.space.genus == 0 and (self.space.n == 4 or self.space.n == 3):            
            return 1
        if self.space.genus == 1 and self.space.n == 1:
            #hard code in these answers
            if self.irr_boundaries == 1:
                return Rational((1,2))
            else: # ch1, psi1, or kappa1, or lambda??
                return Rational((1,24))
        if self.space.genus == 2 and self.space.n == 0:
            if len(self.cherns) == 1 and self.cherns[0] == 3:
                #hard code in this special case
                return Rational((-1, 34560))
                
        return None #meaning it is not a base case.
                
    def key(self):
        """
        Returns a tuple of tuples, suitable for use as a dict key, that represents this monomial up to permutation of marked points.
        
        If there are only psis, returns a simplified tuple compatible with the function from tau.py.
        """
        self.psis.sort(key = lambda t: t[1]) #sort by exponent
        if len(self.red_boundaries) == 0 and len(self.cherns) == 0 and len(self.lambdas) == 0 and self.irr_boundaries == 0 and len(self.kappas) == 0:
            d = from_exp_list_to_tau(self.space.genus, self.space.n, self.psi_exp())            
            return (self.space.genus, tuple(d))
        self.red_boundaries.sort() 
        self.cherns.sort()
        self.lambdas.sort()
        self.kappas.sort()
        reindex_marks = dict( ((psi_tuple[0].mark, i) for i, psi_tuple in enumerate(self.psis) ))
        return ( self.space.genus, self.space.n, tuple( self.psi_exp() ),
                      tuple( self.kappas),
                      tuple( (rb.table_key(reindex_marks) for rb in self.red_boundaries) ),
                      self.irr_boundaries,
                      tuple( self.cherns), 
                      tuple( self.lambdas)
        )
        
    def can_push_down(self):
        """
        Returns true if we can use Faber's trick to push down self-intersections of irreducible boundaries.
        """
        return len(self.red_boundaries) == 0 and len(self.cherns) == 0 and len(self.lambdas) == 0
        
    def push_down_psis(self):
        r"""
        Assume that there are only psis and irrs.  
        
        Returns a polynomial on Mbar_g,0 that is the push down under the forgetful map of the monomial represented by self.
        """
        M = self.space.pushed_down()
        return push_down(self.psi_exp(), self.space.genus, self.space.n) * prod([kappa_class(M,d) for d in self.kappas]) * irreducible_boundary(M)**self.irr_boundaries
        
    def push_down_psis_one(self, mark):
        """
        Compute the intersection number by first pushing down one step using the string equation or the formula (1.7) in [AC96].
        
        Assume that there are only psis and irrs and kappas, and that we have already compensated for pulling back the kappas (i.e. that that kappa_classes in self are really pullbacks of kappa_classes.
        
        This returns a number, since it calls the intersect method itself.
        TODO:  Why?  Does that make sense.        
        """
        M = self.space.forget(mark)
        try:
            index = [psituple[0].mark == mark for psituple in self.psis].index(True) #find the psi class whose mark is being forgotten 
        except ValueError:
            index = None
        
        if index != None:            
            #this is because the mark we want to forget has a psi class, so the string equation does not apply.  We must use (1.7) of [AC96].
            new_kappas = self.kappas + [self.psis[index][1]-1]
            new_psis = list(self.psis)
            new_psis.remove(self.psis[index])
            m = prod([kappa_class(M, d) for d in new_kappas])*prod( [psi_class(M, psi_tuple[0].mark)**psi_tuple[1] for psi_tuple in new_psis]) * irreducible_boundary(M) ** self.irr_boundaries
            return intersect([M], m)
        else:  #the string equation applies (see (1.9) of [AC96])
            value = 0
            for psi_tuple in self.psis:
                value += intersect_monomial_not_product([M],self.monomial.decrease_return_new(psi_tuple[0]))  #optimization potential here!
                #notice that the psis will then not be on the right space... does it matter???
            return value
    

        
    def get_MgnLb_indexes(self):
        r"""
        Returns a list of indexes corresponding to Faber's "MgnLb.txt" program.  This is useful for automated testing purposes.
        """
        l = []
        for tclass, expon in self.monomial.decompose_monomial():
            l += [self.space.num(tclass)]*expon
        return l
        
    def get_a_point(self):
        """
        Returns a marked point that does not have a corresponding psi class if there is such, or any mark otherwise.
        """
        pick_from = list(self.space.marks.difference(set([psituple[1] for psituple in self.psis])))
        if len(pick_from) > 0:
            return pick_from[0]
        else:
            return self.psis[0][1]
         
        
        
class monomial_data_from_product(monomial_data):
    """
    This class is the same as monomial_data except that you have to manually call recieve.  This is useful when dealing with intersections on products, see intersect_monomial_product.
    """
    def __init__(self, space):
        """
        This does not parse the data of the monomial, you have call recieve manually.
        """
        monomial_data.__init__(self, space)
        
    def recieve(self, tclass, expon):
        """
        Notice that this version of recieve also updates self.monomial.
        """
        self.monomial *= tclass**expon
        monomial_data.recieve(self, tclass, expon)

#@breadth_first_tree     
def intersect(spaces, p, check_degree_for_user = False):
    """
    This is the main entry point.
    
    INPUT: 
     - spaces -- A list of Mgn objects that represents the spaces involved, i.e. we are computing an intersection on the product of all the spaces in this list.
     - p -- a polynomial in tautalogical ring classes (see the file TautRing3.py)
     - check_degree_for_user -- if this is true, the function should raise an exception instead of returning 0 if the polynomial contains a monomial of the wrong degree.
     
    OUTPUT:
     - The intersection number for the polynomial. 
    """    
    #print p, spaces#, p.expr
    result = 0
    if len(spaces) == 1:
        intersect_method = intersect_monomial_not_product
    else:
        intersect_method = intersect_monomial_product
    try:
        for m, coeff in p.monomials():
            result += coeff * intersect_method(spaces, m, check_degree_for_user)
        return result
    except AttributeError:
        if isinstance(p, int) or p in QQ:
            tot_dim = sum([space.dimension for space in spaces])
            if tot_dim == 0:
                return p
            else:
                return 0
        else:
            raise  

#@breadth_first_tree  
def intersect_monomial_product(spaces, p, check_degree = False):
    """
    Computed the top intersection number of a monomial on a product space.
    
    INPUT:
     - spaces -- A list of spaces involved in this intersection, i.e. we are computing the intersection on the product of the spaces in this list.
     - p -- A monomial (not a polynomial, and with no coefficient) to compute the intergral of.
     - check_degree -- If true, it should raise an error instead of returning 0 for the wrong degree.
    """
    datas = dict( [(space, monomial_data_from_product(space)) for space in spaces])
    for tclass, expon in p.decompose_monomial():
        datas[tclass.space].recieve(tclass, expon)
    
    ans = 1
    for space in spaces:
        if datas[space].vanishes_by_degree():
            return 0
            
    for space in spaces:
        base_case = datas[space].base_case()
        if base_case != None:
            ans *= base_case
        else:
            ans *= intersect_monomial_with_data(datas[space])
    
    return ans
    
class BadDegreeException(Exception):
    """
    This is the exception to be raised for a wrong degree.
    """
    pass
    
    
#@breadth_first_tree  
def intersect_monomial_not_product(space, m, check_degree = False):
    """
    Computes the intersection of a monomial on a space (not a product space).
    
    INPUT:
     - space -- A length 1 list (to be compatible with the signature of intersect_monomial_product) containing the space to compute the intersection on.
     - m -- A monomial to compute the integral of.
     - check_degree -- Raise an error instead of returning 0 for bad degree.
    """
    data = monomial_data(space[0], m)

    if data.vanishes_by_degree():
        if check_degree:
            raise BadDegreeException("The monomial {0} has degree {1}, while the space {2} has dimension {3}.".format(m, data.degree(), repr(space[0]), space[0].dimension))
        else:
            return 0
    
    base_case = data.base_case()
    if base_case != None:
        return base_case
    else:   
        return intersect_monomial_with_data(data)
    

@remember_convert_args(lambda data: data.key(), master_table)
#@breadth_first_tree#_MgnLb
def intersect_monomial_with_data(data):  
    """
    Computes the itegral based on data from a monomial_with_data object.  
    
    Does not check vanishing or base cases!!! (These should have already been checked)
    
    This is the main implementation of the algorithm.
    """ 
    if len(data.lambdas) != 0:
        to_convert, num_lamb = data.get_lambda()
        data.monomial.decrease_by(to_convert, num_lamb)
        return intersect([data.space], data.monomial * (to_convert.as_chern_chars())**num_lamb)

    if len(data.red_boundaries) != 0:
    
        p1, p2 = data.space.next_marks()
        
        pullback_by = data.get_red() 
        data.monomial.decrease(pullback_by)                   
        
        M1 = pullback_by.component1.add_marks(p1)
        M2 = pullback_by.component2.add_marks(p2)
        
        return Rational((1, pullback_by.degree_of_map)) * intersect( [M1,M2], prod([cl.pullback_red((M1,p1), (M2,p2))**expon for cl, expon in data.monomial.decompose_monomial()]) )
        
    if len(data.cherns) != 0:
            
        chern = data.get_chern()
        data.monomial.decrease(chern)
        value = 0
        # #######
        # Here we implement formula (13) of [Yan10].  Notice that we save a little work by taking advantage of symmety and multiply by 1/degree_of_map instead of 1/2.
        # #######
        for i in range(0, chern.degree): 
            ivalue = 0
            #do the red boundaries
            for psp1, psp2, psi_pair, degree_of_map in chern.get_red_boundaries_and_psis(i): 
                ivalue += Rational((1, degree_of_map)) * intersect([psp1[0], psp2[0]], psi_pair * prod([cl.pullback_red(psp1,psp2)**expon for cl, expon in data.monomial.decompose_monomial()]) )
            
            #do the irr boundary
            (irr_pullback_space, p1, p2), psi_pair = chern.get_irr_boundary_and_psis(i)

            ivalue += Rational((1,2)) * intersect([irr_pullback_space], psi_pair * prod([cl.pullback_irr(irr_pullback_space, p1, p2)**expon for cl, expon in data.monomial.decompose_monomial()]) )
            
            #alternating
            value += (-1)**i * ivalue
        

        value += intersect([data.space], chern.get_kappa()*data.monomial)
        value *= chern.get_coeff()
        
        return value
        
    if data.irr_boundaries > 1:
        #special relation in Faber's paper
        m = max([data.space.genus, 3*data.space.genus - 3])
        if data.irr_boundaries >= m + 1:
            return 0
        # 
    
        if data.space.not_pushed_down(): #use this to avoid self-intersections of irrs on spaces with marked points, speeds it up a bit.
            if len(data.kappas) != 0:
                new_poly = data.monomial
                p = data.get_a_point()
                for d in set(data.kappas):
                    new_poly = new_poly.subs(kappa_class(data.space, d), kappa_class(data.space,d) + psi_class(data.space,p)**d)
                ans = 0
                for m, coef in new_poly.monomials():
                    mdata = monomial_data(data.space,m)
                    ans += coef * mdata.push_down_psis_one(p)  
                return ans
            else:
                return intersect([data.space.pushed_down()], data.push_down_psis())
                
            
    if data.irr_boundaries != 0: #pullback is irreducible boundary
        p1, p2 = data.space.next_marks()   
    
        pullback_by = irreducible_boundary(data.irr_space)
        data.monomial.decrease(pullback_by)
        
        M = Mgn(data.space.genus - 1, data.space.marks.union([p1,p2]) ) 
        
        return Rational((1,2)) * intersect([M], prod([cl.pullback_irr(M, p1, p2)**expon for cl, expon in data.monomial.decompose_monomial()]))    
        

    elif len(data.kappas) != 0:
        p = data.space.next_mark()
        M = data.space.add_marks(p)
        psi_new = psi_class(M,p)
        return intersect([M], psi_new**(data.kappas[0] + 1) * prod( [psi.change_space(M)**expon for psi, expon in data.psis]) * prod( [kappa_class(M, d) - psi_new**(d) for d in data.kappas[1:]]) )  
  
    else:
        return psi_intersect(data.space.genus, data.space.n, data.psi_exp())

