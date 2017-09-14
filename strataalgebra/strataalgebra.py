
from sage.all import cached_method, CommutativeRings, prod, Rational, factorial, floor, Matrix, ZZ, subsets, var, RR, repr_lincomb, vector
from sage.rings.commutative_algebra import CommutativeAlgebra
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import CommutativeAlgebraElement

from allstrata import StrataPyramid
from stratagraph import StrataGraph, PsiKappaVars, graph_count_automorphisms

from intersect import decorate_with_monomial, genericGHstructures, shared_edges

class StrataAlgebraElement(CommutativeAlgebraElement):
    
    def __init__(self,parent,coef_dict):
        if parent is None:
            raise ValueError("The parent must be provided")

        self.coef_dict = coef_dict
            
        CommutativeAlgebraElement.__init__(self, parent)   
        
    def _repr_(self):
        l = self.coef_dict.items()
        l.sort()
        return repr_lincomb([(self._pretty_basis_name(index[0],index[1]), coef) for index, coef in l])
        
        
        #terms = [self._pretty_coef(coef,codim,index) for (codim,index),coef in l if coef != 0]
        #if len(terms) >0:
        #    return " + ".join(terms)
        #return "0"
            
    def _pretty_basis_name(self, codim, index):
        name = self.parent().get_stratum(codim,index).nice_name()
        if name == None:
            name = "s_{0},{1}".format(codim,index)
        return name

        
    def _add_(self, other):
        new_dict = self.coef_dict.copy()
        for (codim,index), coef in other.coef_dict.items():
            dict_plus_equals(new_dict,(codim,index), coef)
        return StrataAlgebraElement(self.parent(), new_dict)
        
    def _mul_(self,other):
        #print self, "*", other
        new_dict = dict()
        for (codim,index), coef in self.coef_dict.items():
            for (codim1,index1), coef1 in other.coef_dict.items():
                for prod_index, prod_coef in self.parent()._prod((codim, index), (codim1, index1)).items():
                    dict_plus_equals(new_dict,(codim+codim1,prod_index),coef*coef1*prod_coef)
        
        return StrataAlgebraElement(self.parent(),new_dict)

    def __div__(self, other):
        #print "lets try to divide!!!!!!!!"
        #print "self", self
        #print "other", other

        if 1/other in self.parent().base():
            return StrataAlgebraElement(self.parent(), {i:c/other for i,c in self.coef_dict.items()})
        else:
            raise Exception("Tried to divide by {0}. Such division not possible!".format(other))

    def integrate(self):
        """
        Return the integral of this class, i.e. its degree in the top codimension.
        
        Classes of codimension less that the dimension of the moduli space will integrate to 0.
        
        This uses the FZ relations to perform the integration. It is probably not very efficient. But it provides a nice check of the implementation. Consider using the ``topintersections`` module if you need to compute something quickly. ::
        
            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,1,(1,2))
            sage: (s.psi(1)*s.psi(2)).integrate()
            1/24
            sage: s.kappa(2).integrate()
            1/24
            sage: s.get_stratum(2,0) #just so you can see it.
            [0 1 2 0 0]
            [0 0 0 2 1]
            [0 1 1 0 1]
            sage: s(2,0).integrate()
            1
            sage: s(1,2).integrate()
            0
            sage: (42*s(2,0) + s(1,1) + 48*s.kappa(2)).integrate()
            44    
        
         
        """
        if self.codim() < self.parent().moduli_dim:
            return 0
        result = 0
        ints = self.parent().basis_integrals()
        for (codim,index), coef in self.coef_dict.items():
            if codim == self.parent().moduli_dim:
                result += coef*ints[index]
        return result

    def dict(self):
        """
        Return a dictionary with keys as :class:`StrataGraph` objects and values as the coefficient of that stratum in this element. ::

            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,1,(1,2,3)); s
            Strata algebra with genus 1 and markings (1, 2, 3) over Rational Field
            sage: a = s.psi(1)*s.psi(2) - 7* s.kappa(3); a
            ps1*ps2 - 7*ka3
            sage: a.dict()
            {ps1*ps2: 1, ka3: -7}

        """
        return {self.parent().strataP.get_stratum(r,i) : c for (r,i), c in self.coef_dict.items() if c !=0}

    def __eq__(self, other):
        if 0 in self.coef_dict.values():
            self.coef_dict = {i:c for i,c in self.coef_dict.items() if c != 0}
        if 0 in other.coef_dict.values():
            other.coef_dict = {i:c for i,c in self.coef_dict.items() if c != 0}
        return self.coef_dict == other.coef_dict
        
    def in_kernel(self):
        """
        Determine whether this :class:`StrataAlgebraElement` is in the span of the FZ relations, and hence in the kernel of the map to the tautological ring. ::
        
            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,0,(1,2,3,4,5))
            sage: b = s.boundary(0,(1,2,5)) + s.boundary(0,(1,2)) - s.boundary(0,(1,3,5)) - s.boundary(0,(1,3))
            sage: b
            Dg0m1_2_5 + Dg0m1_2 - Dg0m1_3_5 - Dg0m1_3
            sage: b.in_kernel()
            True
            sage: (s.psi(1) - s.psi(2)).in_kernel()
            False

            
        It should work fine for non-homogeneous things as well. ::
            
            sage: (b + s.psi(1)).in_kernel()
            False
            sage: (b + s.psi(1)**2 - s.psi(2)**2).in_kernel()
            True
        """
        for codim in range(self.codim()+1):
            v = vector([0]*self.parent().hilbert(codim))
            for (cd, index), coef in self.coef_dict.items():
                if cd == codim:
                    v[index] = coef
            if v not in self.parent().FZ_matrix(codim).row_space():
                return False
        return True
                
            
    def codim(self):
        """
        Returns the codimensions of this :class:`StrataAlgebraElement`. 
        
        If it is not homogeneous, it retunrs the maximum codimension of a basis element with a non-zero coefficient. ::
        
            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,2)
            sage: s(1,0).codim()
            1
            sage: s(2,1).codim()
            2
            sage: s(0,0).codim()
            0
            sage: (35*s(2,3) - 7*s(1,1) + s.kappa(2) - s(0,0)).codim()
            2
            sage: s.zero().codim()
            -1
        """
        codim = -1
        for (cd, index), coef in self.coef_dict.items():
            if cd > codim and coef != 0:
                codim = cd
        return codim
        
        
            
            

        

def dict_plus_equals(d, key, value):
    if key in d.keys():
        d[key]+=value
    else:
        d[key]=value
        
         
     
        
    
class StrataAlgebra(CommutativeAlgebra, UniqueRepresentation):
    """
    StrataAlgebra
    """

    #Need these to make the autodocs work.
    get_stratum = StrataPyramid.get_stratum
    FZ_matrix = StrataPyramid.FZ_matrix
    FZ_matrix_pushforward_basis = StrataPyramid.FZ_matrix_pushforward_basis
    print_strata = StrataPyramid.print_strata

    def __init__(self,base,g,markings=(), make_vars = True):
        r"""
        A ring representing the Strata algebra.

        :param Ring base: The ring of coefficients you want to work over, usually ``QQ``.
        :param int g: The genus.
        :param tuple markings: The markings should be positive integers. Repeats are allowed. Defaults to no markings.
        :param bool make_vars: Defaults to True. If True, creates variables ``ps``, ``ps_``, ``ka1``, ... , ``ka{d}`` (where d is the dimension of the moduli space) that can be used to create basis elements.
            
        First import the module: ::
        
            sage: from strataalgebra import *

        Construct a :class:`StrataAlgebra`: ::
            
            sage: SA = StrataAlgebra(QQ,1,(1,2)); SA
            Strata algebra with genus 1 and markings (1, 2) over Rational Field

        Print the basis elements in a certain codimension, with their (arbitrary) index: ::

            sage: SA.print_strata(2)
            **** i: 0
            [0 1 2 0 0]
            [0 0 0 2 1]
            [0 1 1 0 1]
            <BLANKLINE>
            **** i: 1
            [0 1 2 0 0]
            [0 0 1 1 1]
            [0 1 0 1 1]
            <BLANKLINE>
            **** i: 2
            [  0   1   2   0]
            [ka1   1   1   2]
            <BLANKLINE>
            **** i: 3
            [     0      1      2      0]
            [     0      1      1 ps + 2]
            <BLANKLINE>
            **** i: 4
            [     0      1      2      0]
            [     0 ps + 1      1      2]
            <BLANKLINE>
            **** i: 5
            [     0      1      2      0]
            [     0      1 ps + 1      2]
            <BLANKLINE>
            **** i: 6
            [      0       1       2       0]
            [      0       1       1       1]
            [ka1 + 1       0       0       1]
            <BLANKLINE>
            **** i: 7
            [     0      1      2      0]
            [     0      1      1      1]
            [     1      0      0 ps + 1]
            <BLANKLINE>
            **** i: 8
            ka2
            <BLANKLINE>
            **** i: 9
            ka1^2
            <BLANKLINE>
            **** i: 10
            ka1*ps1
            <BLANKLINE>
            **** i: 11
            ka1*ps2
            <BLANKLINE>
            **** i: 12
            ps1^2
            <BLANKLINE>
            **** i: 13
            ps1*ps2
            <BLANKLINE>
            **** i: 14
            ps2^2
            <BLANKLINE>


        Classes that are monomials in :math:`\psi` and :math:`\kappa` have self-explanatory names.

        More complicated classes are displayed in a matrix as follows:

        Each row after the first corresponds to a vertex of the graph.
        The constant term in the entry in the first column is the genus.
        The kappa classes also appear in the first column.
        Each column beyond the first corresponds to an edge or a half edge.
        The entry in the first row gives the label of the half edge, or ``0`` for a full edge.
        The constant term of the entry in location (v,e) gives the number of times (0, 1, or 2) that edge e touches vertex v.
        A ``ps`` in entry (v,e) means a :math:`\psi`-class associated to the half edge coming out of v.
        A ``ps_`` may occur when there is a loop at the vertex.
        The entry in the top left is just padding.

        To create classes, you can use their codimension and index.
        Boundary strata and psi-kappa monomials are represented by their special names,
        and the rest are represented by ``s_{codim},{index}`` ::

            sage: a = SA(2, 1); a
            s_2,1
            sage: b = SA(2, 7); b
            s_2,7
            sage: c = SA(2, 11); c
            ka1*ps2
            sage: d = SA(1,1); d
            Dg1

        Vector space arithmetic is supported. ::

            sage: 3*b
            3*s_2,7
            sage: a*72 - b/17
            72*s_2,1 - 1/17*s_2,7

        Use :meth:`~strataalgebra.StrataAlgebra.get_stratum` if you need to know what an unamed basis element means. ::

            sage: SA.get_stratum(2,7)
            [     0      1      2      0]
            [     0      1      1      1]
            [     1      0      0 ps + 1]

        You can construct :math:`\psi,\;\kappa` monomials and boundary divisors with the methods :meth:`~strataalgebra.StrataAlgebra.kappa`, :meth:`~strataalgebra.StrataAlgebra.psi`,
        :meth:`~strataalgebra.StrataAlgebra.boundary`, and :meth:`~strataalgebra.StrataAlgebra.irr`.

        You can construct an element using the matrix notation. Just pass a list of lists into your :class:`StrataAlgebra`. ::
        
            sage: var('ps') #This is usually done automatically, but we have to do it manually here for the dotests to work.
            ps
            sage: SA([[0,1,2,0],[1,0,0,ps+1],[0,1,1,1]])
            s_2,7

        Here is an example of the ``ps_``: ::

            sage: s = StrataAlgebra(QQ,2,())
            sage: s.get_stratum(3,9)
            [       0        0]
            [       1 ps^2 + 2]
            sage: s.get_stratum(3,10)
            [         0          0]
            [         1 ps*ps_ + 2]
            sage: var('ps_')
            ps_
            sage: s([[0,0],[1,ps_*ps+2]])
            s_3,10

        One of the main features is the computation of the product. ::

            sage: a*b
            0

        Of course, the codimension was too big. Lets do some less trivial ones. ::

            sage: SA = StrataAlgebra(QQ,1,(1,2,3)); SA
            Strata algebra with genus 1 and markings (1, 2, 3) over Rational Field
            sage: SA.psi(1) * SA.psi(2)
            ps1*ps2
            sage: SA.get_stratum(2,7) #just so you can see what it is
            [0 1 2 3 0 0]
            [0 0 1 1 1 0]
            [0 1 0 0 1 1]
            [1 0 0 0 0 1]
            sage: SA(2,7)*SA.psi(1)
            0
            sage: SA(2,7)*SA.kappa(1)
            s_3,36
            sage: SA.get_stratum(3,36)
            [      0       1       2       3       0       0]
            [      0       0       1       1       1       0]
            [      0       1       0       0       1       1]
            [ka1 + 1       0       0       0       0       1]
            sage: SA.irr()^3
            6*s_3,4 - 3*s_3,15 - 3*s_3,23 - 3*s_3,35 + s_3,48 + s_3,49

        Everything should work with distributive laws, etc. ::

            sage: SA.kappa(2)*(3*SA.psi(3)-SA(1,3)) == -SA(1,3) *SA.kappa(2) + SA.psi(3)*SA.kappa(2)*3
            True

        It should work over any ring. ::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: SA = StrataAlgebra(R, 2); SA
            Strata algebra with genus 2 and markings () over Univariate Polynomial Ring in t over Integer Ring
            sage: (3+t)*SA.kappa(1) + t^2 + t*SA.kappa(1)
            t^2*one + (2*t+3)*ka1
            sage: _^2
            t^4*one + (4*t^3+6*t^2)*ka1 + (4*t^2+12*t+9)*ka1^2


        There may be problems over a non-divisible ring. ::

            sage: SA.irr()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Rational Field' and 'Strata algebra with genus 2 and markings () over Univariate Polynomial Ring in t over Integer Ring'

        Also, repeated names of the markings are allowed. The following corresponds to :math:`\overline{\mathcal M}_{1,2}/S_2`.
        Compare the codimension 2 strata to our earlier example. ::

            sage: SA = StrataAlgebra(QQ,1,(1,1)); SA
            Strata algebra with genus 1 and markings (1, 1) over Rational Field
            sage: SA.print_strata(2)
            **** i: 0
            [0 1 1 0 0]
            [0 0 0 2 1]
            [0 1 1 0 1]
            <BLANKLINE>
            **** i: 1
            [0 1 1 0 0]
            [0 0 1 1 1]
            [0 1 0 1 1]
            <BLANKLINE>
            **** i: 2
            [  0   1   1   0]
            [ka1   1   1   2]
            <BLANKLINE>
            **** i: 3
            [     0      1      1      0]
            [     0      1      1 ps + 2]
            <BLANKLINE>
            **** i: 4
            [     0      1      1      0]
            [     0 ps + 1      1      2]
            <BLANKLINE>
            **** i: 5
            [      0       1       1       0]
            [      0       1       1       1]
            [ka1 + 1       0       0       1]
            <BLANKLINE>
            **** i: 6
            [     0      1      1      0]
            [     0      1      1      1]
            [     1      0      0 ps + 1]
            <BLANKLINE>
            **** i: 7
            ka2
            <BLANKLINE>
            **** i: 8
            ka1^2
            <BLANKLINE>
            **** i: 9
            ka1*ps1
            <BLANKLINE>
            **** i: 10
            ps1^2
            <BLANKLINE>
            **** i: 11
            ps1*ps1
            <BLANKLINE>
        """
    
        if base not in CommutativeRings():
            raise ValueError("%s is no CommutativeRing" % base)
            
        self.g = g 
        self.markings = markings
        
        self.moduli_dim = 3*g-3+len(markings)
            
        self.strataP = StrataPyramid(g,markings)
                
        self._prod_dict = dict()
        
        self.get_stratum = self.strataP.get_stratum
        self.FZ_matrix = self.strataP.FZ_matrix
        self.print_strata = self.strataP.print_strata
        self.FZ_matrix_pushforward_basis = self.strataP.FZ_matrix_pushforward_basis
        
        CommutativeAlgebra.__init__(self, base)

        if make_vars:
            if self.moduli_dim == 1:
                var("ka1")
            elif self.moduli_dim > 1:
                var(["ka{0}".format(a) for a in range(1, self.moduli_dim+1)])
            var(StrataGraph.ps_name,StrataGraph.ps2_name)

        
    def _repr_(self):
        return "Strata algebra with genus {0} and markings {1} over {2}".format(self.g, self.markings, self.base())
    
    
    
    @cached_method #This breaks the autodocs for some reason...
    def basis_integrals(self, s=None):
        """
        Return a list of numbers corresponding to the integrals of the basis elements in the top codimension.

        This is computed via the FZ relations, so it is probably not fast. However, it is a nice check.

        The value is cached, so you only have to compute it once per session.

        This is used by :meth:`~strataalgebra.StrataAlgebraElement.integrate` which is the more likely way you will want to use it.
        ::

            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,1,(1,))
            sage: s.print_strata(1)
            **** i: 0
            D_irr
            <BLANKLINE>
            **** i: 1
            ka1
            <BLANKLINE>
            **** i: 2
            ps1
            <BLANKLINE>
            sage: s.basis_integrals()
            [1, 1/24, 1/24]

        """
        int_kappa = Rational((1, 24**self.g * factorial(self.g)))
        G = StrataGraph(self.strataP.open_stratum().M)
        G.M[1, 0] += StrataGraph.Xvar ** self.moduli_dim

        kappa_index = self.strataP.index_from_graph_codim(G,self.moduli_dim) #._dstrata[self.moduli_dim][G]

        M = self.FZ_matrix(self.moduli_dim)
        ker = M.right_kernel()
        if ker.dimension() != 1:
            raise Exception("Something unexpected here!")
        return list( int_kappa / ker.gen(0)[kappa_index] * ker.gen(0))
        
            
    def base_ring(self):
        return self.base().base_ring()
        
    def characteristic(self):
        return self.base().characteristic()
        
    def do_all_products(self):
        """
        Compute all the products of basis elements of the ring (thus caching their values).
        This could take a long time for some rings. ::
        
            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,1,(1,))
            sage: s.do_all_products()

        """
        for r1 in range(1, self.moduli_dim):
            for r2 in range(1, self.moduli_dim-r1+1):
                for i in range(self.hilbert(r1)):
                    for j in range(self.hilbert(r2)):
                        print "******** doing", (r1,i),(r2,j)
                        self._prod((r1, i), (r2, j))
        
    def _prod(self, ci1, ci2):
        """
        Compute a product. You should probably not call this directly. Instead, just use the * operator. Values are cached in ``self._prod_dict``.

        INPUT:

        - ``ci1`` -- A tuple of (codim,index), where codim is the codimension of the basis element and index is its index in the list of classes.
        - ``ci2`` -- Ditto.

        OUTPUT: A dictionary mapping indexes of basis elements in codimension ci1[1]+ci2[1] to coefficients.
        """
        #print ci1, ci2
        if ci1[0] + ci2[0] > self.moduli_dim:
            #print "codim wrong!!"
            return dict()
            
        if ci1 > ci2:
            temp = ci1
            ci1 = ci2
            ci2 = temp

        if ci1 == (0,0):
            #print "identity!"
            return {ci2[1]: 1}
        
        value = self._prod_dict.get((ci1, ci2), None)
        if value != None:
            #print "retrived value", value
            return value
                
        #print "computing!"
        
        #compute the value!
        G = self.strataP.get_stratum(*ci1) #[ci1[0]][ci1[1]]
        H = self.strataP.get_stratum(*ci2) #[ci2[0]][ci2[1]]
       # print "got stratum"
        
        prod_dict = dict()
        loop_coef = 2**(G.num_undecorated_loops() + H.num_undecorated_loops())
        
        for A in self.strataP.common_degenerations(G,H):
            PsiKappaRing, kappa, psi, psi2 = PsiKappaVars(A)
            
            autA = graph_count_automorphisms(A)
            
            
            for Gs,Hs in genericGHstructures(G,H,A):
                #print "A:"
                #print A
                #print "Gs:"
                #print G 
                #print Gs
                #print "Hs:"
                #print H
                #print Hs
                
                #print "psi", psi
                #print "psi2", psi2
                
                prod_f_A_G = self._prod_f_A_G(Gs, kappa, psi, psi2)
                prod_f_A_H = self._prod_f_A_G(Hs, kappa, psi, psi2)
                
                excess = 1
                for v1,v2,eind in shared_edges(Gs,Hs):
                    #print "psi", psi
                    #print "psi2", psi2
                    if v1 != v2:
                        excess *=(-psi[(eind,v1)]-psi[(eind,v2)])
                    else:
                        excess *=(-psi[(eind,v1)]-psi2[(eind,v2)])
                #print "did excess"
                        
                F_A_GH  = prod_f_A_G*prod_f_A_H*excess
                #print prod_f_A_G, "*" , prod_f_A_H, "*", excess, "=", F_A_GH
                #print "F_A_GH", F_A_GH
                #print "A"
                #print A
                #print
                
                #print "aut", autA
                #print "loopc", loop_coef
                
                prod_codim = ci1[0] + ci2[0]
                #print prod_codim
                
                
                
                if F_A_GH in ZZ:
                    if F_A_GH == 0:
                        continue
                    #print autA, loop_coef, F_A_GH/autA*loop_coef
                    dict_plus_equals(prod_dict, self.strataP.index_from_graph_codim(A,prod_codim), Rational((F_A_GH * loop_coef, autA)))
                else:
                    for m,c in F_A_GH.dict().items():
                        Ad = decorate_with_monomial(PsiKappaRing,A,m)
                        Ad_index = self.strataP.index_from_graph_codim(Ad,prod_codim,None) #_dstrata[prod_codim].get(Ad, None)
                        #print "Ad:"
                        #print Ad
                        #print "Ad_index", Ad_index, c/autA*loop_coef
                        if not Ad_index is None:
                            #print "!!!!!!!!!!!!!"
                            #print autA
                            dict_plus_equals(prod_dict, Ad_index, Rational((c * loop_coef,autA)))
                    
        self._prod_dict[(ci1,ci2)] = prod_dict
        return prod_dict

    @staticmethod
    def _prod_f_A_G(Gs, kappa, psi, psi2):
        """
        Used in the computation of product. Do not call directly.
        """
        result = 1
        G = Gs.G.strataG
        A = Gs.A.strataG
        for v in range(1,G.num_vertices()+1):
            #print "Gs", Gs
            #print "G,g",v
            result *= prod(( sum((kappa.get((j,w),0) for w in Gs.alpha_inv[v]))**f for j,f in G.kappa_on_v(v))) * prod(( Gs.psi_no_loop(edge,v,ex,psi) for edge, ex in G.psi_no_loop_on_v(v) )) * prod(( Gs.psi_loop(edge,v, ex1, ex2, psi,psi2) for edge, ex1, ex2 in G.psi_loop_on_v(v) ))

        return result                    
        
    
        
        
    def _element_constructor_(self, *args, **kwds):
        """
        See __init__ documentation.
        """

        if len(args) == 1:
            if args[0] in self.base():
                return StrataAlgebraElement(self,{(0,0):args[0]})
            elif isinstance(args[0], StrataGraph):
                G = args[0]
            else:
                G = StrataGraph.from_nice_matrix(args[0])
            return StrataAlgebraElement(self, {self.strataP.codim_index_from_graph(G):1})

        if len(args) == 2:
            return StrataAlgebraElement(self,{(args[0],args[1]):1})
        
    def _coerce_map_from_(self,S):
        if self.base().has_coerce_map_from(S):
            return True
        if isinstance(S,StrataAlgebra):
            if self.base().has_coerce_map_from(S.base()):
                return True 
        return False
    
    def psi(self,mark):
        r"""
        Return a psi class.
        
        :param int mark:  The mark that the :math:`\psi`-class is associated to.
        :rtype: :class:`StrataAlgebraElement`

        Examples ::
        
            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,0,(1,2,3,4,5))
            sage: psi1 = s.psi(1); psi1
            ps1
            sage: psi2 = s.psi(2); psi2
            ps2
            sage: psi1*psi2
            ps1*ps2
            sage: (psi1*psi2).integrate()
            2

        Note that the psi variables are not automatically injected into the namespace ::

            sage: ps1
            Traceback (most recent call last):
            ...
            NameError: name 'ps1' is not defined

        In case of repeated marks, notice that the :math:`\psi`-class is the pushforward under the quotient map. Observe: ::

            sage: s = StrataAlgebra(QQ,1,(1,1))
            sage: s.psi(1)^2
            ps1^2 + ps1*ps1
            sage: s.psi(2)
            Traceback (most recent call last):
            ...
            ValueError: tuple.index(x): x not in tuple
            sage: s = StrataAlgebra(QQ,1,(1,1,1))
            sage: s.psi(1)^2
            2*ps1^2 + 4*ps1*ps1
            sage: s.psi(1)^3
            4*ps1^3 + 24*ps1^2*ps1 + 8*ps1*ps1*ps1
            sage: var('ps') #This should be done automatically.
            ps
            sage: s([[0,1,1,1],[1,ps^3+1,1,1]])
            ps1^3
            sage: s = StrataAlgebra(QQ,0,(1,1,1,2,2))
            sage: s.psi(1)*s.psi(2)
            12*ps1*ps2
            sage: s.psi(1)*s.psi(1)
            4*ps1^2 + 8*ps1*ps1

        This maybe looks surprising, but it makes sense with the formula

        .. MATH ::

            \pi_*(\alpha) \pi_*(\beta) = \frac{1}{|G|} \pi_* \left( \sum_{\sigma, \tau \in G} \sigma_*(\alpha) \tau_*(\beta) \right)

        for two class :math:`\alpha,\;\beta \in H^*(X)` where :math:`\pi: X \rightarrow X/G` is the quotient map and :math:`G` is a finite group.
        """
        #G = StrataGraph(self.strataP._dstrata[0].keys()[0].M)
        G = StrataGraph(self.strataP.open_stratum().M)
        index = self.markings.index(mark)+1
        G.M[1,index] += StrataGraph.Xvar
        G.compute_invariant()
        return self(G)
        
    def kappa(self, a):
        """
        Return a kappa class.

        :param int a: The subscript (codimension) of the kappa class you want.
        :rtype: :class:`StrataAlgebraElement`

        ::

            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,2,())
            sage: s.kappa(1)
            ka1
            sage: s.kappa(3)
            ka3
            sage: s.kappa(1)*s.kappa(2)
            ka1*ka2
            sage: s.kappa(2)*s.kappa(2)
            0
            sage: s.boundary(1,())*s.kappa(1)
            s_2,4
            sage: s.get_stratum(2,4)
            [      0       0]
            [ka1 + 1       1]
            [      1       1]


        The subscript must be less than or equal to the dimension of the moduli space. ::

            sage: s.kappa(4)
            Traceback (most recent call last):
            ...
            KeyError: 4
        """
        #G = StrataGraph(self.strataP._dstrata[0].keys()[0].M)
        G = StrataGraph(self.strataP.open_stratum().M)
        G.M[1,0] += StrataGraph.Xvar**a
        #G.compute_invarient()
        return self(G)
        
    def boundary(self, g1, markings1 = ()):
        """
        Return a boundary divisor with genus `g1` and `markings1` points on one component

        :param int g1: The genus on one component
        :param list markings1: A list or tuple of markings on the one component. Defaults to the empty list.
        :rtype: :class:`StrataAlgebraElement`

        Note that if the two components are the same, it will return the stratum with a coefficient of 1/2. ::
        
        
            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ, 1, (1,2,3))
            sage: s.boundary(0, (1,2))
            Dg0m1_2
            sage: s.boundary(1,[3])
            Dg0m1_2
            sage: s = StrataAlgebra(QQ, 2)
            sage: s.boundary(1)
            1/2*Dg1

        Each component must still be stable. ::

            sage: s.boundary(2)
            Traceback (most recent call last):
            ...
            KeyError: Dg0

        """
        markings2 = list(self.markings)
        for m in markings1:
            markings2.remove(m)
        g2 = self.g-g1
        M = Matrix(StrataGraph.R, [[-1, 0] + list(markings1) + markings2, [g1,1] + [1]*len(markings1) + [0]*len(markings2), [g2,1] + [0]*len(markings1) + [1]*len(markings2)])
        G = StrataGraph(M)
        G.compute_invariant()
        
        if self.g -g1 == g1 and len(self.markings) == 0: 
            return Rational((1,2))* self(G)
        
        return self(G)
        
    def irr(self):
        """
        Make the irreducible boundary. It will be returned with a coefficient of 1/2. ::
        
            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ, 1, (1,2))
            sage: s.irr()
            1/2*D_irr
        """
        # type: () -> object
        if self.g < 1:
            raise Exception("No irreducible boundaries for genus 0!")
        #G = StrataGraph(self.strataP._dstrata[0].keys()[0].M)
        G = StrataGraph(self.strataP.open_stratum().M)
        G.add_edge(1,1)
        G.M[1,0]-=1
        G.compute_invariant()
        return Rational((1,2))*self(G)
        
    def MgnLb_int(self, index_list):
        """
        Computes an integral of boundary divisors, kappa classes, and psi classes.

        :param index_list: A list of indices of classes, according the the scheme of Carl Faber's ``MgnLb`` Maple program. This function is useful because so you can test our implementation of the product and the FZ_relations.
        :rtype: `Rational`

        Examples: ::

            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,1,(1,2))
            sage: s.MgnLb_int([1,6])
            1/2
            sage: s.MgnLb_int([1,2])
            1/24
            sage: s.MgnLb_int([4])
            1/24

        .. SEEALSO ::

            :meth:`~strataalgebra.StrataAlgebra.MgnLb_class`

        """
        return prod(( self.MgnLb_class(i) for i in index_list )).integrate()
        
    def MgnLb_class(self,index):
        """
        Returns the class corresponding to the index from Carl Faber's ``MgnLb`` Maple program.
        This is useful for testing purposes. ::

            sage: from strataalgebra import *
            sage: s = StrataAlgebra(QQ,1,(1,2))
            sage: s.MgnLb_class(1)
            ps1
            sage: s.MgnLb_class(4)
            ka2
            sage: s.MgnLb_class(6)
            1/2*D_irr
            sage: s.MgnLb_class(2)
            ps2

        """
        #print "making classes again!"
        if index <= len(self.markings):
            return self.psi(index)
        index -= len(self.markings)
        if index <= self.moduli_dim:
            return self.kappa(index)
        index -= self.moduli_dim
        if index <= self.g:
            raise Exception("We don't do the ch classes!")
        index -= self.g
        if index == 1:
            return self.irr()
        index -=1
        
        marks = set(self.markings)
        reducible_boundaries = []
        if len(self.markings) != 0:
            first_mark_list = [marks.pop()]            
            for g1 in range(0, self.g + 1):
                for p in subsets(marks):
                    r_marks = set(first_mark_list + p)
                    if 3*g1 - 3 + len(r_marks) + 1 >= 0 and 3*(self.g-g1) - 3 + len(self.markings) - len(r_marks) + 1 >= 0:
                        reducible_boundaries.append( (g1, r_marks) )  
                        
            reducible_boundaries.sort(key = lambda b: sorted(list(b[1])))
            reducible_boundaries.sort(key = lambda b: len(b[1]))
            reducible_boundaries.sort(key = lambda b: b[0])
        
        else: #self.n == 0
            for g1 in range(1, floor(self.g/2.0)+1):
                reducible_boundaries.append( (g1, [])) 
            
        return self.boundary(*reducible_boundaries[index-1])   
        
    def hilbert(self, codim = None):
        """
        Give the number of basis elements in the Strata algebra for the given codimension (the Hilbert function).
        If the codimension is omitted, return a list of all of them.

        This is NOT the Betti numbers for the moduli space of curves! For that see :meth:`~strataalgebra.StrataAlgebra.FZ_betti`.

        :param int codim: Optional. The codimension you want.

        ::

            sage: from strataalgebra import *
            sage: StrataAlgebra(QQ, 2, (1,)).hilbert()
            [1, 4, 17, 49, 92]
            sage: StrataAlgebra(QQ, 2, (1,)).hilbert(2)
            17


        """
        if codim is None:
            return [len(self.strataP.all_strata(i)) for i in range(self.moduli_dim + 1)]
        return len(self.strataP.all_strata(codim))

    def FZ_betti(self, codim = None):
        """
        Give the dimension of the cohomology of moduli space of curves, as predicted by the Faber-Zagier relations.
        If the codimension is omited, return a list of all of them.

        :param codim: Optional. The codimension you want.

        ::
            
            sage: from strataalgebra import *
            sage: StrataAlgebra(QQ, 2, (1,)).FZ_betti()
            [1, 3, 5, 3, 1]
            sage: StrataAlgebra(QQ, 2, (1,)).FZ_betti(2)
            5

        .. SEEALSO ::
            :meth:`~strataalgebra.StrataAlgebra.hilbert`

        """
        if codim is None:
            return [self.FZ_betti(i) for i in range(self.moduli_dim+1)]
        else:
            return len(self.strataP.all_strata(codim)) - self.strataP.FZ_matrix_pushforward_basis(codim).rank()



#yH = StrataGraph(
#Matrix(StrataGraph.R,[[-1,0,0,0,0,0],[2,1+X,0,0,0,0],[0,1,1,1,0,0],[0,0,1,0,2,0],[0,0,0,1,0,2]]))
#yH.compute_invariant()


#yG = StrataGraph(Matrix(StrataGraph.R,[[-1,0],[3+X^2,2]]))
#yG.compute_invariant()

