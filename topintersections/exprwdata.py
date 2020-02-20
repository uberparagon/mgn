"""
This is my work-around for not being able to assign arbitrary attributes to elements of the sage Symbolic Ring, or to inherit from symbolic ring element.

We create a class ExprWithDataGen that represents a variable.  The class ExprWithData stored a pair of dictionaries that created a bijection between ExprWithDataGen instances and Symbolic Ring variables that it creates for this purpose.  An instance of the class ExprWithData represents an expression in ExprWithDataGens.

ExprWithData also has some methods that are useful for breaking down polynomials into their terms and monomials into factors.
"""
from __future__ import print_function

from sage.all import *


class ExprWithDataGen(object):
    """
    Reprents an atomic expression in your ring of things with data.  This is the subclass of TautRingElement.
    """
    def __add__(self, other):
        return other + ExprWithData.get_gen(self)
        
    __radd__ = __add__
    
    def __mul__(self, other):
        return other * ExprWithData.get_gen(self) 
        
    __rmul__ = __mul__
    
    def __pow__(self, other):
        return ExprWithData.get_gen(self)**other
        
    def __neg__(self):
        return ExprWithData.get_gen(self)*(-1)
        
    def __sub__(self, other):
        return self + (-1)*other
        
    def __div__(self, other):
        return self * (Rational(1)/other)
        
    def monomials(self):
        yield ExprWithData.get_gen(self), 1
        
    @property
    def expr(self):
        """
        Returns the underlying expression.
        """
        return ExprWithData.get_gen(self).expr
        


class ExprWithData(object):
    """
    Represents an expression involving combinations of ExprWithDataGens.
    """
    
    var_obj_dict = dict()
    obj_var_dict = dict()
    var_name_prefix = "VAR_"
    var_name_suffix = "_RAV"
    
    @classmethod
    def reset(cls):
        """
        Resets the dictionaries, in case they get too big.  Not sure if this is really necessary, but you could try it if memory seems to be a problem.  
        
        Of course any existing instances of ExprWithData will become unusable after calling this.
        """
        cls.var_obj_dict = dict()
        cls.obj_var_dict = dict()
        #lambda_class.as_chern_chars_dict = dict()
    
    def __init__(self, expr):
        """
        INPUT:
         - expr -- A SymbolicRing expression.  The symbols in the expression should have already been in the dictionaries.  You should probably use arithmetic operators if you want to build an expression, don't call this directly.
        """
        self.expr = expr
        
    @classmethod 
    def get_gen(cls, gen):
        """
        Pass in an ExprWithDataGen (e.g. a psi_class, a kappa_class) and get a variable from the symbolic ring that will be linked to it using the dictionaries.  If a corresponding variable already exists, this method just returns it.
        """
        v = cls.obj_var_dict.get(gen) 
        if  v != None:
            return ExprWithData(v)
            
        var_ = SR.var(cls.var_name_prefix + str(len(cls.var_obj_dict)) + cls.var_name_suffix)        
        cls.var_obj_dict[var_] = gen
        cls.obj_var_dict[gen] = var_
        
        return ExprWithData(var_)
        
    def is_constant(self):
        """
        NOTE:: this is surprizingly slow, if I recall.  Avoid this in inner loops if possible.
        """
        return self.expr in CC
        

    def __repr__(self):
        """
        This is important... turns the arbitrary variable names into the names you really want.
        """
        repr_str = repr(self.expr)
        for v, o in self.var_obj_dict.items():
            repr_str = repr_str.replace(repr(v), repr(o))
        return repr_str
            
    def __add__(self, other):
        if isinstance(other, ExprWithDataGen):
            return other + self
        if isinstance(other, ExprWithData):
            return ExprWithData(self.expr + other.expr)
        #otherwise, assume other is a constant
        return ExprWithData(self.expr + other)
        
    def __neg__(self):
        return ExprWithData(-self.expr)
        
    def __sub__(self, other):
        return self + (-1)*other
        
    __radd__ = __add__
        
    def __mul__(self, other):
        if isinstance(other, ExprWithDataGen):
            return other * self
        if isinstance(other, ExprWithData):
            return ExprWithData(self.expr * other.expr)
        #assume other in QQ or something
        return ExprWithData(self.expr * other)
        
    def __div__(self,other):
        return self * Rational((1,other))
        
    __rmul__ = __mul__
    
    def __pow__(self, other):
        return ExprWithData(self.expr**other)
        
    def expand(self):
        print("Warning, monomials already expands it!")
        return ExprWithData(self.expr.expand())
        
    def monomials(self):
        """
        Gives a generator.

        Yields tuples of (monomial, coefficient) for this expression.
        """
        expr = self.expr.expand()
        if expr.operator() == operator.pow or expr.operator() ==None:
            if expr.is_zero():
                return
            yield self, 1
            return
        
        if expr.operator() == sage.symbolic.operators.mul_vararg:
            gothrough = [expr]
        elif expr.operator() == sage.symbolic.operators.add_vararg:
            gothrough = expr.operands()
        
        for m in gothrough:
            if m.operator() == None:
                if m in QQ:
                    yield ExprWithData(SR(1)), m
                else:
                    yield ExprWithData(m), 1
            elif m.operator() == operator.pow:
                yield ExprWithData(m),1
            else:                    
                if m.operands()[-1] in QQ:
                    coeff = m.operands()[-1]
                else:
                    coeff = 1
                yield ExprWithData(m/coeff), coeff
                
                
    def decrease(self, item):
        """
        Reduces the exponent of the given item by 1.
        """
        self.expr = self.expr/self.obj_var_dict[item]
        
    def decrease_by(self, item, n):
        """
        Reduces the exponent of the given item by n.
        """
        self.expr = self.expr/(self.obj_var_dict[item])**n
    
    def decrease_return_new(self, item):
        """
        Creates a new ExprWithData, reduces the exponent of the given item, and returns the new one.
        """
        return ExprWithData(self.expr/self.obj_var_dict[item])
        
    def decompose_monomial_check_1(self):
        """
        Checks if this expression is 1, and if so, returns an empty list.  Otherwise, just calls decompose_monomial.
        
        NOTE::  The check self.expr == 1 is somewhat expenseive, and often not neccesary, which explains why have this method instead of putting the check in decompose_monomial.
        """
        if self.expr == 1:
            return []
        else:
            return self.decompose_monomial()
        
        
    def decompose_monomial(self):
        """
        Generator. Yeilds pairs (class, exponent).  Assumes that the expression is a monomial.  For example, psi1^2*kappa1^3 will yeild
        (psi1,2)
        (kappa1,3) 
        """
        if self.expr.operator() == operator.pow or self.expr.operator() == None:
            gothrough = [self.expr]
        elif self.expr.operator() == sage.symbolic.operators.mul_vararg:
            gothrough = self.expr.operands()
        else:
            raise Exception("Not a monomial!")
              
            
        for v in gothrough:
            if v.operator() == operator.pow:
                expon = v.operands()[1]
                o = self.var_obj_dict[v.operands()[0]]
            else:
                expon = 1
                o = self.var_obj_dict[v]
            yield o, int(expon)
            
    def subs(self, sub_out, sub_in):
        """
        Does a substitution, see implementation.
        """
        return ExprWithData(self.expr.subs(sub_out.expr == sub_in.expr))
            

        
#Some stuff to test this out if you want to play with it.

# class myVar(ExprWithDataGen):
    # def __init__(self, name):
        # self.name = name
        
    # def __repr__(self):
        # return self.name
        
    # def __hash__(self):
        # return hash(self.name)
        
    # def __eq__(self, other):
        # return self.name == other.name
        
# a = myVar("a")
# b = myVar("b")
# c = myVar("c")
# d = myVar("d")
# e = myVar("e")
