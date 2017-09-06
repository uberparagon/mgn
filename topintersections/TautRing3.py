try:
    from sage.all import *
    from exprwdata import *
    from LazyAttr import lazy_attr
    from richcomparisonmixin import RichComparisonMixin
except ImportError:
    pass

class Mgn(object):
    """
    Represents the compactification of the moduli space of curves with genus g and specified marked points.
    
    Sometimes this is also used to represent the data associated to a reducible boundary divisor.  From one point of view, I think I regret that design decision, but I don't think it will affect performace or maintainability too much.
    """
    def __init__(self, genus, marks, var_scope = None):
        """
        INPUT:
         - genus
         - marks -- Either an integer, or a list or set of the labels for the marks
         - var_scope -- A place to create the variables for classes on this space.  Mostly you would pass in globals(), I think.  If you omitt this argument, it doesn't create the variables.
        """
        try:
            self.marks = frozenset(marks)
        except TypeError:
            self.marks = frozenset(range(1, marks+1))
        self.genus = genus
        
        if var_scope != None:
            self.mak(var_scope)
            
        self.hash_value = hash(self.marks) + hash(self.genus)
        #if 3*self.genus - 3 + self.n < 0:
        #    raise Exception, "No such space exists! g= " +str(genus) + " marks = " + str(marks)
        
    def __eq__(self, other):
        return self.marks == other.marks and self.genus == other.genus
        
    def __hash__(self):
        return self.hash_value #hash(self.marks) + hash(self.genus)
        
    @property
    def dimension(self):
        """Gives the dimension of the space."""
        return 3*self.genus - 3 + self.n
        
    def add_marks(self, *marks):
        """Returns a new space that is this one plus some extra marks."""
        return Mgn(self.genus, self.marks.union(marks))
        
    @property
    def n(self):
        """Returns the number of marked points in this space"""
        return len(self.marks)
        
    def num(self, cls):
        """Returns Faber's index for this class, if it is one of the MgnLb ones.  Useful for automated testing."""
        try:
            return self.classes.index(cls) + 1
        except:
            raise Exception, "Someone asked for " + str(cls) + ", which is not here " + str(self)
            
        
    def complementary_component(self, big_space):
        """
        Assuming that this Mgn represents a reducible componet, then find its complement in the big_space.        
        """
        if not self.marks.issubset(big_space.marks):
            print self, big_space
            raise Exception, "Bad marked points"
        if not self.genus <= big_space.genus:
            raise Exception, "Bad genus"
        return Mgn(big_space.genus - self.genus, big_space.marks.difference(self.marks))
        
        
    def same_as(self, space, in_space):
        """
        Determines whether self and space represent the same reducible boundary in in_space. Returns true or false.
        """
        if self.marks == space.marks and self.genus == space.genus:
            return True
        space = space.complementary_component(in_space)
        if self.marks == space.marks and self.genus == space.genus:
            return True
        return False
        
    def __repr__(self):
        return "Mbar_" + str(self.genus) + "_" + str(self.n)
        
    def __str__(self):
        return "Compactified moduli space of genus {0} curves with {1} marked points {2}.".format(self.genus, self.n, str(self.marks)[11:-2])
        
        
    @lazy_attr
    def classes(self):
        """
        Returns a list of Tautalogical classes on this space.  The order should correspond to that of Faber's program.
        """
        #print "making classes again!"
        l = []
        for p in self.marks:
            l.append(psi_class(self,p))
        for d in range(1, self.dimension + 1):
            l.append(kappa_class(self,d))
        for i in range(1, self.genus+1):
            l.append(chern_char(self, 2*i-1))
        if True:#self.genus != 0:
            l.append(irreducible_boundary(self))
        marks = set(self.marks)
        reducible_boundaries = []
        if self.n != 0:
            first_mark_list = [marks.pop()]            
            for g1 in range(0, self.genus + 1):
                for p in subsets(marks):
                    r_marks = set(first_mark_list + p)
                    if 3*g1 - 3 + len(r_marks) + 1 >= 0 and 3*(self.genus-g1) - 3 + self.n - len(r_marks) + 1 >= 0:
                        reducible_boundaries.append( reducible_boundary(self, Mgn(g1, r_marks))  )
                        
            reducible_boundaries.sort(key = lambda b: sorted(list(b.component1.marks)))
            reducible_boundaries.sort(key = lambda b: len(b.component1.marks))
            reducible_boundaries.sort(key = lambda b: b.component1.genus)
        
        else: #self.n == 0
            for g1 in range(1, floor(self.genus/2.0)+1):
                reducible_boundaries.append(reducible_boundary(self, Mgn(g1, []))) 
                
        
        l += reducible_boundaries            
        
        for i in range(1,self.genus+1):
            l.append(lambda_class(self,i))
        return l
        
    def __getitem__(self, index):
        """
        Accessor method to get the classes.  Shifts the index so it matches Fabers.  This index is printed by the rij() function.
        """
        if index in self.marks:
            return psi_class(self, index)
        return self.classes[index-1]
    
    def rij(self):        
        """
        Prints a list of all the Tautalogical classes for this space, with their index.  The name is borrowed from Faber, I think it is Dutch.  I provide an engligh alias, but I don't know a short name to call it.
        """
        for i, c in zip(range(len(self.classes)), self.classes):
            print "[" + str(i+1) + "]  " + repr(c) #+ " --- " + str(c) 
            
    print_classes = rij
            
    def mak(self, where):
        """
        Creates symbols in the namespace ``where`` for all the classes on this space.  This makes it so you can type in e.g. 
        
            sage: psi1*ka2
        
        and it knows what you mean.
        """
        for long_name_for_class_that_probably_wont_be_in_global_namespace in self.classes:
            if long_name_for_class_that_probably_wont_be_in_global_namespace != 0:
                exec(repr(long_name_for_class_that_probably_wont_be_in_global_namespace) + " = long_name_for_class_that_probably_wont_be_in_global_namespace", locals(), where)
            
    def red_boundaries_as_spaces(self):
        """
        Generator object, gives all the reducible boundaries for this space, but returns them as spaces, i.e., as a Mgn object with the defining genus and marked points.
        
        This is used to compute itersections with ch.
        """
        marks = set(self.marks)
        if self.n != 0:
            first_mark_list = [marks.pop()]
       
            
            p1,p2 = self.next_marks()     
                
            for g1 in range(0, self.genus + 1):
                for p in subsets(marks):
                    r_marks = set(first_mark_list + p)
                    if 3*g1 - 3 + len(r_marks) + 1 >= 0 and 3*(self.genus-g1) - 3 + self.n - len(r_marks) + 1 >= 0:
                        yield  (Mgn(g1, r_marks.union([p1])), p1),   (Mgn(self.genus - g1, marks.difference(r_marks).union([p2])), p2) 
                    
        else: # self.n == 0
            for g1 in range(1, floor(self.genus/2.0)+1):
                yield (Mgn(g1, [1]), 1) , (Mgn(self.genus-g1, [2]), 2)
                
                    
    def irr_boundary_as_space(self):
        """
        Gives the irreducible boundaries for this space, but returns it as a space, i.e., as a Mgn object with the defining genus and marked points.
        
        This is used to compute itersections with ch.
        """
        p1,p2 = self.next_marks()
        
        return Mgn(self.genus - 1, self.marks.union([p1,p2])), p1, p2 

    def next_marks(self):
        """
        Returns the two next marks that haven't been used on this space.
        """
        if self.n != 0:
            pmax = max(self.marks)
        else:
            pmax = 0
            
        p1 = pmax + 1
        p2 = pmax + 2
        
        return p1, p2
        
    def next_mark(self):
        """
        Returns the next mark that hasn't been used.
        """
        if self.n != 0:
            pmax = max(self.marks)
        else:
            pmax = 0
            
        return pmax + 1
        
        
    def __iter__(self):
        """
        To catch a common mistake...
        """
        raise Exception, "Don't iterate this!  Did you pass this to intersect without putting it in a list?"
        
        
    def degree_index_dict(self):
        """
        Returns a dictionary, with keys being the degrees, and values being a list of indexes of classes of the corrsponding degree.  This is used by a different function to generate a list of all potentially non-zero monomials.
        
        Does not return psi classes, since they are made differently.
        """
        did = dict()
        for i,c in enumerate(self.classes):
            if isinstance(c, lambda_class) or isinstance(c, psi_class) or c == 0:
                continue            
            try:
                degree = c.degree
            except AttributeError:
                degree = 1
            if not did.has_key(degree):
                did[degree] = []
            did[degree].append(i+1)
        return did
        
    def degree_index_dict_no_red(self):
        """
        Returns a dictionary, with keys being the degrees, and values being a list of indexes of classes of the corrsponding degree.  This is used by a different function to generate a list of all potentially non-zero monomials.  This version does do the lambdas, but no reducible boundaries or psis.
        """
        did = dict()
        for i,c in enumerate(self.classes):
            if isinstance(c, reducible_boundary) or isinstance(c, psi_class)  or c == 0:
                continue            
            try:
                degree = c.degree
            except AttributeError:
                degree = 1
            if not did.has_key(degree):
                did[degree] = []
            did[degree].append(i+1)
        return did
        
    def not_pushed_down(self):
        """
        Returns true if there is a forgetful map out of this space.
        """
        return (self.genus >= 2 and self.n != 0) or (self.genus == 1 and self.n > 1) or (self.genus == 0 and self.n > 3)
        
    def pushed_down(self):
        """
        Returns the space obtained by forgeting as many marked points as possible.
        """
        if self.genus == 0:
            return Mgn(0,3)
        if self.genus == 1:
            return Mgn(1,1)
        return Mgn(self.genus,0)
        
    def forget(self, mark):
        """
        Returns a new space with the specified point forgotten.
        """
        return Mgn(self.genus, self.marks.difference([mark]))
  

        
class TautRingElement(ExprWithDataGen):
    """
    Base class of all Tautalogical ring classes.
    """
    
    def __init__(self, space): 
        ExprWithDataGen.__init__(self)
        self.space = space
        
    def __new__(cls, space, *args):
        if space.dimension < 0:
            #print "Warning, attempted to create class on bad space {0}, returning 0 instead!".format(space)
            return 0
        return ExprWithDataGen.__new__(cls)
    
        
    def pullback_red(self, (M1, p1), (M2, p2)):
        """
        This method must be overridden.
        
        Pull this class back to the product M1 x M2, with p1 and p2 being the new marks.
        """
        pass
        
    def pullback_irr(self, M, p1, p2):
        """
        This method must be overridden.
        
        Pull back to the irreducible curve, where the class is currently on M, and p1 and p2 are the names of the new marks.
        """
        pass
        
        


class psi_class(TautRingElement):
    def __init__(self, space, mark):
        TautRingElement.__init__(self,space)
        self.mark = mark   
        
    def change_space(self, space):
        return psi_class(space, self.mark)
        
    def pullback_red(self, (M1, p1), (M2, p2)):
        if self.mark in M1.marks:
            return self.change_space(M1)
        else:
            return self.change_space(M2)
            
    def pullback_irr(self, M, p1, p2):
        return self.change_space(M)
        
    def __repr__(self):
        return "psi" + str(self.mark)
        
    def __str__(self):
        return ("Psi class for marked point {0} on " + repr(self.space)).format(self.mark)
        
    def __eq__(self, other):
        return isinstance(other, psi_class) and self.space == other.space and self.mark == other.mark
        
    def __hash__(self):
        return hash(self.space) + hash(self.mark)
        
class nice_class(TautRingElement):
    """
    Base class for kappa_class, lambda_class, and chern_char, since they have similar properties.
    """
    def __init__(self, space, degree):
        TautRingElement.__init__(self,space)
        self.degree = degree
        
    def __new__(cls, space, degree):
        if degree > space.dimension:
            return 0
        else:
            return TautRingElement.__new__(cls, space)
        
    def change_space(self, new_space):
        return self.__class__(new_space, self.degree)
        
    def pullback_red(self, (M1, p1), (M2, p2)):
        #NOTE: this is not correct for higher lambdas!!!!!!!!!!
        return self.change_space(M1) + self.change_space(M2)
        
    def pullback_irr(self, M, p1, p2):
        return self.change_space(M)
        
    def __repr__(self):
        return self.name + str(self.degree)
        
    def __str__(self):
        return "{0} class of degree {1} on {2}".format(self.long_name.capitalize(), self.degree, repr(self.space))
        
    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.space == other.space and self.degree == other.degree
    
    def __hash__(self):
        return hash(self.space) + hash(self.degree) + hash(self.name)
        
class kappa_class(nice_class):
    name = "ka"
    long_name = "kappa"
    def __new__(self, space, degree):
        if degree == 0:
            return space.genus*2 - 2 + space.n
        return nice_class.__new__(kappa_class, space, degree)
        

        
class reducible_boundary(TautRingElement, RichComparisonMixin):

    def __init__(self, space, component1):
        TautRingElement.__init__(self,space)
        if self.space.n != 0:
            if sorted(self.space.marks)[0] in component1.marks:
                self.component1 = component1
            else:
                self.component1 = component1.complementary_component(space)
        else:
            if 2*component1.genus > space.genus:
                self.component1 = component1.complementary_component(space)
            else:
                self.component1 = component1
            

        
    def __new__(cls, space, component1):
        if component1.genus > space.genus:
            #print "Warning, tried to create a bad (genus too big) reducible boundary on {0} represented by {1}".format(space, component1)
            return 0
        if component1.dimension < -1:
            #print "Warning, tried to create a bad reducible boundary on {0} represented by {1}".format(space, component1)
            return 0
        if component1.complementary_component(space).dimension < -1:
            #print "Warning, tried to create a bad reducible boundary (complement is bad) on {0} represented by {1}".format(space, component1)
            return 0
        return TautRingElement.__new__(cls, space)
        
        
    def __eq__(self, other):
        return isinstance(other, reducible_boundary) and self.space == other.space and self.component1 == other.component1
        
    def __hash__(self):
        return hash(self.space) + hash(self.component1)
        
    @property
    def component2(self):
        return self.component1.complementary_component(self.space)
        
    def table_key(self, reindex_dict):
        """
        Returns a tuple that can be used to represent this class in the dictionary.  It recieves as a argument a dictionary that tells you how to reindex the marks.
        """
        reindexed_marks = []
        for m in self.component1.marks:
            new_m = reindex_dict.get(m)
            if new_m == None:
                if len(reindex_dict) == 0:
                    new_m = 0
                else:
                    new_m = max(reindex_dict.values())+1
                reindex_dict[m] = new_m
            reindexed_marks.append(new_m)
        return tuple( [self.component1.genus] + sorted(reindexed_marks) )
        
    def change_space(self,space):
        #print "change space", self, space
        return reducible_boundary(space, self.component1)
        
    def pullback_red(self, (M1, p1), (M2, p2)):
        """
        See [Fab99].
        """
        if self.component1.add_marks(p1) == M1: 
            return_value = - psi_class(M1, p1) - psi_class(M2, p2)
            if self.space.n == 0 and M1.genus < M2.genus:
                return return_value + reducible_boundary(M2, Mgn(M2.genus - M1.genus, [p2]))
            elif self.space.marks.issubset(M1.marks) and M1.genus >= M2.genus > 0:
                return return_value + reducible_boundary(M1, Mgn(M1.genus - M2.genus, M1.marks))
            else:
                return return_value              
        
        #else, not a self intersection:
        if self.space.n == 0:
            if M1.genus > self.component1.genus:
                return reducible_boundary(M1, Mgn(M1.genus - self.component1.genus, [p1])) + reducible_boundary(M2, Mgn(M2.genus - self.component1.genus, [p2]))
            else:
                if 2*self.component1.genus == self.space.genus:
                    return reducible_boundary(M2, Mgn(M2.genus - self.component1.genus, [p2]))
                else:
                    return reducible_boundary(M2, Mgn(M2.genus - self.component1.genus, [p2])) + reducible_boundary(M2, Mgn(self.component1.genus - M1.genus, [p2]))
               
               
        #else if n != 0:        
        return_value = 0
        if M1.genus >= self.component1.genus and self.component1.marks.issubset(M1.marks):
            return_value += reducible_boundary(M1, self.component1)
        if self.component1.genus >= M1.genus and M1.marks.difference([p1]).issubset(self.component1.marks):
            return_value += reducible_boundary(M2, Mgn(self.component1.genus - M1.genus, self.component1.marks.difference(M1.marks).union([p2])))
        if self.component1.genus >= M2.genus and M2.marks.difference([p2]).issubset(self.component1.marks):
            return_value += reducible_boundary(M1, Mgn(self.component1.genus - M2.genus, self.component1.marks.difference(M2.marks).union([p1])))
            
        return return_value
        

        
                   
    
    def pullback_irr(self, M, p1, p2):
        """
        See [Fab99].
        """
        if len(self.space.marks) == 0 and 2 * self.component1.genus == self.space.genus:
                return self.change_space(M)
        else:
            return reducible_boundary(M, Mgn(self.component1.genus - 1, self.component1.marks.union([p1,p2]))) + self.change_space(M) 
                        
    @property
    def degree_of_map(self):
        """
        Returns either 1 or 2, depending on the degree of the map from the product space to this reducible boundary.  It should be two if there is an automorphism of the graph, and 1 otherwise.
        """
        if 2*self.component1.genus == self.space.genus and self.space.n == 0:
            return 2
        else:
            return 1
        
    def __repr__(self):
        if self.component1.n != 0:
            return "Dg{0}m{1}".format(self.component1.genus, str(self.component1.marks)[11:-2].replace(", ", "_"))
        else:
            return "Dg{0}".format(self.component1.genus)
        
    def __str__(self):
        return ("Boundary class on " + repr(self.space) + " corresponding to genus {0} and points {1}.").format(self.component1.genus, str(self.component1.marks)[11:-2])
        
    #These methods are for sorting, i.e. to decide which to pull back by first.  Not sure if this is optimal, but heuristically seems approximately correct.
    def comp1(self):
        return abs(2*self.component1.genus-self.space.genus) + abs(2*self.component1.n-self.space.n)
    def comp2(self):
        return abs(2*self.component1.n-self.space.n)
    def comp3(self):
        return self.component1.genus
        
        
    def __lt__(self,other):
        for cmp in (reducible_boundary.comp1, reducible_boundary.comp2, reducible_boundary.comp3):
            if cmp(self) < cmp(other):
                return True
            if cmp(self) > cmp(other):
                return False
        return sorted(list(self.component1.marks)) < sorted(list(other.component1.marks))
        
        
class irreducible_boundary(TautRingElement):

    def __new__(cls, space):
        if space.genus == 0:
            #print "Trying to create a irreducible boundary on {0}, returning 0 instead.".format(space)
            return 0
        return TautRingElement.__new__(cls, space)

    def change_space(self, space):
        return irreducible_boundary(space)
        
    def pullback_red(self, (M1, p1), (M2, p2)):
        """
        See [Fab99].
        """
        return irreducible_boundary(M1) + irreducible_boundary(M2)
        
    def pullback_irr(self, M, p1, p2):  
        """
        See [Fab99].
        """
        return_value=  -psi_class(M, p1) - psi_class(M, p2) + irreducible_boundary(M)
        if self.space.n == 0:
            for h in range(1, self.space.genus -2 + 1):
                return_value += reducible_boundary(M, Mgn(h, [p1]))
            return return_value
        else:
            for h in range(0, self.space.genus - 1 +1):
                nbar = set(self.space.marks)
                point1 = nbar.pop() #now nbar is missing an element
                for Marks_minus1 in subsets(nbar):
                    if set(Marks_minus1) == nbar and h == self.space.genus - 1:
                        continue
                    Marks = set(Marks_minus1 +[point1])

                    return_value += reducible_boundary(M, Mgn(h, Marks.union([p1]))) + reducible_boundary(M, Mgn(h, Marks.union([p2])))
            return return_value 
            

        
    def __eq__(self, other):
        return isinstance(other, irreducible_boundary) and self.space == other.space
        
    def __hash__(self):
        return hash(self.space)
        
    def __repr__(self):
        return "irr"
        
    def __str__(self):        
        return ("Irreducible boundary class on " + repr(self.space))
        
class chern_char(nice_class):
    name = "ch"
    long_name = "chern character"
        
    def __new__(cls, space, degree):
        if space.genus == 0 or degree > 2*space.genus-1:
            return 0
        else:
            return TautRingElement.__new__(cls,space)
        
    def get_coeff(self):
        """
        Returns the coefficent in formula that coverts chern classes to kappas and psis (formula for chern character in [Yan10]).
        """
        return  bernoulli(self.degree+1) / factorial(self.degree + 1) 
        
    def get_kappa(self):
        """
        Returns the "bad" (Mumford's original) kappa of the same degree.
        """
        return kappa_class(self.space, self.degree) - sum((psi_class(self.space, p)**self.degree for p in self.space.marks) )
        
    def get_red_boundaries_and_psis(self,i):
        """
        See formula for chern character in [Yan10].  This is a generator that yeilds a tuple of data used to compute the sum over the reducible boundary pushforward maps.
        """
        for (M1,p1), (M2, p2) in self.space.red_boundaries_as_spaces():
            if M1.n == 1 and M2.n == 1 and M1.genus == M2.genus:
                degree_of_map = 2
            else:
                degree_of_map = 1
            yield (M1, p1), (M2, p2), psi_class(M1, p1) ** i * psi_class(M2, p2) ** (self.degree - i - 1), degree_of_map
            
    def get_irr_boundary_and_psis(self, i):
        """
        See formula for chern character in [Yan10].  This returns a tuple of data used to compute the sum over the irreducible boundary pushforward map.
        """
        M, p1, p2 = self.space.irr_boundary_as_space()
        return (M, p1, p2), psi_class(M, p1) ** i * psi_class(M, p2) ** (self.degree - i - 1) 

class lambda_class(nice_class):
    name = "la"
    long_name = "lambda"
    def __new__(cls, space, degree):
        if degree > space.genus:
            return 0
        return nice_class.__new__(cls, space, degree)
    
    #as_chern_chars_dict = dict()
    
    def pullback_red(self, (M1, p1), (M2, p2)):
        raise Exception("Pullback has not been implemented for lambdas, but you should be able to avoid this.")
        
    def as_chern_chars(self):
        """
        Returns an equivilant expression terms of chern characters.
        """
        #value = self.as_chern_chars_dict.get(self)
        if True: #value == None:
            
            #print "on space", self.space, "with lambda_", self.degree
            upto = ceil(self.degree/2.0)+1
            
            S = PowerSeriesRing(SR, "t")
            f = ( sum( [factorial(2*i - 2) * ExprWithData.get_gen(chern_char(self.space, 2*i-1)).expr * S.gen()**(2*i-1) for i in range(1, upto)] ) ).exp(self.degree+1)
            #print f
            #print f.coefficients()
            value = ExprWithData(f[self.degree]) #wow, way faster to use indexing rather than coefficients list
            #self.as_chern_chars_dict[self] = value
        return value       
