import collections
import itertools
from stratagraph import StrataGraph, graph_count_automorphisms, aut, graph_isomorphic, PsiKappaVars, X_to_kappa
from sage.all import Rational, Integer, cached_method, cached_function, IntegerVectors, Partitions, Matrix, Permutations, binomial, factorial, floor, zero_matrix
from intersect import decorate_with_monomial
from sage.structure.unique_representation import UniqueRepresentation

A_list = [factorial(6*n)/(factorial(3*n)*factorial(2*n)) for n in range(100)]
B_list = [factorial(6*n+1)/((6*n-1)*factorial(3*n)*factorial(2*n)) for n in range(100)]


class OrderedSetWithIndex(collections.OrderedDict):
    
    #def __init__(self):
    #    super(OrderedSetWithIndex,self).__init__()
    #    self.cur_index = 0
        
        
    def add(self, item):
        if item in self.keys():
            return
        collections.OrderedDict.__setitem__(self, item, len(self.keys()))
        
    def add_all(self, iter):
        for i in iter:
            self.add(i)
     
class StrataPyramidMetaclass(type):
    def __init__(self, *args, **kwargs):
        print "Depreciated!!!!"
        super(StrataPyramidMetaclass, self).__init__(*args, **kwargs)
        self.dict = dict()
        
    def __call__(self, *args, **kwargs):
        g = args[0]
        markings = tuple(args[1])
        
        value = self.dict.get((g,markings),None)
        if value == None:
            print "Building StrataPyramid for g={0} and markings {1}".format(g,markings)
            sp = super(StrataPyramidMetaclass,self).__call__(g,markings)
            sp.build()
            self.dict[(g,markings)] = sp
            return sp
        else:
            return value
      
        
class StrataPyramid(UniqueRepresentation):
    #__metaclass__ = StrataPyramidMetaclass
    def __init__(self,genus,markings):
        self.moduli_dim = 3*genus-3 + len(markings)
        if self.moduli_dim < 0:
            raise Exception("Wrong dim!")
        self.genus = genus
        self.markings = markings
        self._specialization = dict() #maps codimension to a set of undecorated strata.
        self._dstrata = dict()
        #self._specialization = dict() #maps a graph to the set of its specializations.
        for r in range(self.moduli_dim+1):
            self._specialization[r] = None #collections.OrderedDict()
            self._dstrata[r] = None #OrderedSetWithIndex()
        #self.FZ_coeff_dict = dict()

        #build level 0 in the constructor.
        G = StrataGraph()
        G.add_vertex(self.genus)
        for m in self.markings:
            G.add_edge(1, 0, m)
        #G.compute_invariant()
        self._specialization[0] = collections.OrderedDict()
        self._specialization[0][G] = set()
        self._dstrata[0] = OrderedSetWithIndex()
        self._dstrata[0].add(G)


    def open_stratum(self):
        return self._dstrata[0].keys()[0]

    def build(self):
        raise Exception("Depreciated!")
        self.build_undecorated_strata()
        self.build_decorated_strata()
          
    def build_undecorated_strata(self):
        raise Exception("Depreciated!")
        
        for r in range(self.moduli_dim):
            for G in self._specialization[r]:
                Gd = degenerate1_list(G)
                self._specialization[r][G] = set(Gd)
                self._specialization[r+1].update(( (G,None) for G in Gd ))
        
        for topG in self._specialization[self.moduli_dim].keys():
            self._specialization[self.moduli_dim][topG] = set()
        print "undecorated construction complete!"

    def codim_index_from_graph(self, G):
        codim = G.codim()
        if self._dstrata[codim] is None:
            self.build_decorations(codim)
        return (codim, self._dstrata[codim][G])

    def index_from_graph_codim(self,G,codim, default = "error"):
        if self._dstrata[codim] is None:
            self.build_decorations(codim)
        value = self._dstrata[codim].get(G,None)
        if value is None:
            if default == "error":
                raise Exception("Graph not found!")
            else:
                return default
        else:
            return value

    def build_specialization(self, codim):
        for r in range(1, codim+1):
            if not self._specialization[r] is None:  # already built
                continue
            #print self.genus, self.markings,"Building specialization level",r
            self._specialization[r] = collections.OrderedDict()

            for G in self._specialization[r-1]:
                Gd = degenerate1_list(G)
                self._specialization[r-1][G] = set(Gd)
                if r == self.moduli_dim:
                    self._specialization[r].update(( (G, set()) for G in Gd))
                else:
                    self._specialization[r].update(( (G, None) for G in Gd))


    def build_decorations(self, codim):
        if self._specialization[codim] is None:
            self.build_specialization(codim)
        for r in range(codim + 1):
            if not self._dstrata[r] is None:
                continue

            #print self.genus, self.markings, "decorating strata codim", r
            self._dstrata[r] = OrderedSetWithIndex()
            self._dstrata[r].add_all(self._specialization[r].keys())
            for less in range(1,r+1):
                for G in self._specialization[r-less]:
                    self._dstrata[r].add_all(decorate1_list(G, less))

    def print_ddims(self):
        self.build_decorations(self.moduli_dim)
        for r in range(self.moduli_dim+1):
            print r, len(self._dstrata[r].keys())
      
    def get_stratum(self, r, j):
        """
        Get the :class:`StrataGraph` associated a a codimension and an index.

        :param r: The codimension
        :param j: The index
        :rtype: :class:`StrataGraph`

        See :class:`StrataAlgebra` documentation for examples.
        """
        if self._dstrata[r] is None:
            self.build_decorations(r)
        return self._dstrata[r].keys()[j]
        
    def print_strata(self,r):
        """
        Print all the strata, with their indexes, in codimension `r`.
        :param r: The codimension
        :return: None

        See :class:`StrataAlgebra` documentation for examples.
        """
        if self._dstrata[r] is None:
            self.build_decorations(r)

        for i,g in enumerate(self._dstrata[r].keys()):
            print "**** i:",i
            print g
            print
    
    @cached_method        
    def _convert_kappa_basis(self,r):
        """
        Returns a matrix that when multiplied on the right of the FZ_matrix with pushforward of psi basis will give the FZ_matrix for the kappa monomial basis.
        """
        M = zero_matrix(len(self._dstrata[r]))
        for Gind,G in enumerate(self._dstrata[r].keys()):
            #print G
            #Gf = forget_decorations(G)
            R, kappa = PsiKappaVars(G,kappa_only = True)
            kappa_dec = 1
            for v in range(1, G.num_vertices()+1):
                #print kappa, G.M[v,0]
                kappa_dec *= X_to_kappa(G.M[v,0],lambda a: kappa[(a,v)])
                #print G.M[v,0], X_to_kappa(G.M[v,0],lambda a: kappa[(a,v)])
            
            #print "kd", kappa_dec
            if kappa_dec == 1:
                M[Gind,Gind] = 1
            else:
                for mon, coef in kappa_dec.dict().items():
                    Gd = decorate_with_monomial(R, G.forget_kappas(), mon)
                    M[Gind, self._dstrata[r][Gd]] += coef
                
        return M #.transpose()
        
    def build_decorated_strata(self):
        """
        You should have called build_decorated_strata first.
        """
        raise Exception("depreciated!")
        #if self._specialization[r] is None:
        #    self.build_specialization(r)

        for r in range(self.moduli_dim+1):
            self._dstrata[r].add_all(self._specialization[r].keys())
            for G in self._specialization[r]:
                for d in range(1, self.moduli_dim-r+1):
                    self._dstrata[r + d].add_all(decorate1_list(G, d))
        print "decorated construction complete!"
                
    def specialization(self,G):
        codim = G.codim_undecorated()
        if codim == self.moduli_dim:
            return set()
        if self._specialization[codim+1] is None:
            self.build_specialization(codim+1)
        return self._specialization[codim][G]

                     
    def all_degenerations_set(self,Gset):
        result = set()
        for Gi in Gset:
            result.update(self.specialization(Gi))
        return result
                                    
    def common_degenerations(self,G,H):
        result = set()
        G = G.forget_decorations()
        H = H.forget_decorations()
        if G.num_edges() > H.num_edges():
            temp = G
            G=H
            H=temp
        
        
        Gd = [G]
        for r in range(H.codim()-G.codim()):
            Gd = self.all_degenerations_set(Gd)
            #print "r=",r, len(Gd)

        if H in Gd:
            result.add(H)
            #Gd.remove(H)
            
        Hd = [H]
        for j in range(0, 2*G.codim()): #check this
            Gd = self.all_degenerations_set(Gd)
            Hd = self.all_degenerations_set(Hd)
            #print "j=",j,len(Gd),len(Hd)
            common = Gd.intersection(Hd)
            result.update(common)
            #Hd.difference_update(common) #you can't do this, I guess
            #Gd.difference_update(common)
            #print "j=",j, "len(Gd)", len(Gd), "len(Hd)", len(Hd)
        
        return result

    def all_strata(self,codim):
        if self._dstrata[codim] is None:
            self.build_decorations(codim)
        return self._dstrata[codim].keys()
    
    def FZ_matrix(self,r):
        """
        Return the matrix of Faber-Zagier-Pixton relations.

        :param r:
        :rtype: :class:`Matrix`

        The columns correspond to the basis elements of the Strata algebra, and each row is a relation.

        Notice that this matrix considers the kappa classes to be in the monomial basis. Thus, is different than the
        output of Pixton's original `taurel.sage` program.

        .. SEEALSO ::

            :meth:`~StrataAlgebra.FZ_matrix_pushforward_basis`


        """
        return Matrix(self.list_all_FZ(r))*self._convert_kappa_basis(r)

    def FZ_matrix_pushforward_basis(self,r):
        """
        Return the matrix of Faber-Zagier-Pixton relations, using the "pushforward" basis, NOT the kappa monomial basis
        that the rest of the code uses.

        :param r:
        :rtype: :class:`Matrix`

        The columns correspond to the basis elements of the Strata algebra, and each row is a relation.
        This matrix should be the same as Pixton's original `tautrel.sage` program after permuting columns.
        """
        return Matrix(self.list_all_FZ(r))

    @cached_method    
    def list_all_FZ(self, r):
        ###
        markings = self.markings
        g=self.genus
        #moduli_type = 4
        ###
        generators = self._dstrata[r].keys() #= get_all_strata(g,r,markings,moduli_type)
        ngen = len(generators)
      
        relations = []
        FZpl = FZ_param_list(3*r-g-1,markings)
        #print "%s codim 0 relations to compute" % len(FZpl)
        #print "FZpl", FZpl
        ccccc = 0
        for FZ_param in FZpl:
      #  if ccccc % 5 == 0:
      #    print "%s done" % ccccc
            ccccc += 1
            relation = [FZ_coeff(G,FZ_param) for G in generators]
            #print "rel", relation
            relations.append(relation)
        old_count = len(relations)

        for r0 in range(1,r):
            strata = self._dstrata[r0].keys()#get_all_strata(g,r0,markings,moduli_type)
            for G in strata:
                vertex_orbits = graph_count_automorphisms(G,True)
                for i in [orbit[0] for orbit in vertex_orbits]:
                    #print "i", i
                    good = True
                    for j in range(G.M.ncols()):
                        if G.M[i,j][0] != G.M[i,j]: #removed an explicit cast here
                            good = False
                            break
                    if good:
                        g2 = G.M[i,0][0]
                        if 3*(r-r0) < g2 + 1:
                            continue
                        d = G.degree(i)
                        if dim_form(g2,d) < r-r0:
                            continue
                        strata2P = StrataPyramid(g2, tuple(range(1,d+1)))
                        #strata2P.build()
                        strata2 = strata2P.all_strata(r-r0) #strata2P._dstrata[r - r0].keys() #get_all_strata(g2,r-r0,tuple(range(1,d+1)),moduli_type) #need help here!!!!!!!
                        #which_gen_list = [-1 for num in range(len(strata2))]
                        which_gen_dict = dict()
                        for G2 in strata2:
                            G_copy = StrataGraph(G.M)
                            G_copy.replace_vertex_with_graph(i,G2)
                            G_copy.compute_invariant()
                            G_copy_index = self._dstrata[r].get(G_copy, None)
                            #print "sanity"
                            #print G_copy
                            #print
                            #print self.get_strata(r,G_copy_index)
                            #print
                            
                            if G_copy_index != None:
                                which_gen_dict[G2] = G_copy_index 
                            #k = self.dstrata[r][G_copy]
                            #for k in range(ngen):
                            #    if graph_isomorphic(G_copy,generators[k]):
                            #        which_gen_list[num] = k
                            #        break
                        rFZpl = reduced_FZ_param_list(G,i,g2,d,3*(r-r0)-g2-1)
                        #print "Computing %s relations to insert at vertex %s into" % (len(rFZpl), i)
                        #print "which"
                        #for G2, num2 in which_gen_dict.items():
                        #    print G2
                        #    print "|-->"
                        #    print self.get_strata(r,num2)
                        #    print
                        #    print
                        ccccc = 0
                        for FZ_param in rFZpl:
                            relation = [0 for k in range(ngen)]
                            #for num in range(len(strata2)):
                            for G2, num2 in which_gen_dict.items():
                                #if which_gen_list[num] != -1:
                                    relation[num2] += FZ_coeff(G2,FZ_param) #FZ_coeff(FZ_param,num,g2,r-r0,tuple(range(1,d+1)),moduli_type)
                            relations.append(relation)
                            ccccc += 1
              #  if ccccc % 5 == 0:
              #    print "%s done" % ccccc
      #print "%s relations in all" % len(relations)
        if len(relations) == 0:
            relations.append([0 for i in generators])
        return relations
    
    
        
    
    
  
@cached_function  
def FZ_coeff(G,FZ_param):
    sigma = FZ_param[0]
    marking_vec = FZ_param[1]
    #g = self.genus
    #markings = self.markings
    
    #print "entering FZ_coeff"
    #print sigma, marking_vec
    #print G
    ####
    #G = self.get_strata(r,num) #get_single_stratum(num,g,r,markings,moduli_type)
    nv = G.num_vertices()
    graph_auts = graph_count_automorphisms(G) #get_autom_count(num,g,r,markings,moduli_type)
    h1_factor = 2**G.h1()
    num_parities = 2**nv
    target_parity = sum((G.parity(i+1) << i) for i in range(nv))

    marking_factors = FZ_marking_factor(G,marking_vec)
    kappa_factors = FZ_kappa_factor_graph(G,sigma) #get_FZ_kappa_factor(num,sigma,g,r,markings,moduli_type)
    hedge_factors = FZ_hedge_factor(G) #get_FZ_hedge_factor(num,g,r,markings,moduli_type)
    #print "mf", marking_factors
    #print "kf", kappa_factors
    #print "hf", hedge_factors
    #print "h1_f", h1_factor
    #print "graph_auts", graph_auts
    #print "np", num_parities
    #print "tp", target_parity
    
    total = Integer(0)
    for i in range(num_parities):
        if marking_factors[i] == 0:
            continue
        for j in range(num_parities):
            total += marking_factors[i]*kappa_factors[j]*hedge_factors[i ^ j ^ target_parity] #changed ^^ to ^

    total /= h1_factor*graph_auts
    
    #print "returning", total
    return total

@cached_function        
def FZ_hedge_factor(G):
  nv = G.num_vertices()
  num_parities = 2**nv
  ne = G.num_edges()
  edge_list = []
  for k in range(1,ne+1):
    if G.M[0,k] == 0:
      edge_list.append([k])
      for i in range(1,nv+1):
        if G.M[i,k] != 0:
          edge_list[-1].append(i)
        if G.M[i,k][0] == 2:
          edge_list[-1].append(i)
  hedge_factors = [0 for i in range(num_parities)]
  for edge_parities in itertools.product(*([0,1] for i in edge_list)):
    parity = 0
    for i in range(len(edge_list)):
      if edge_parities[i] == 1:
        parity ^= 1 << (edge_list[i][1]-1) #supposed to be xor here
        parity ^= 1 << (edge_list[i][2]-1)
    hedge_factor = 1
    for i in range(len(edge_list)):
      if edge_list[i][1] == edge_list[i][2]:
        hedge_factor *= dual_C_coeff(G.M[edge_list[i][1],edge_list[i][0]][1],G.M[edge_list[i][1],edge_list[i][0]][2],edge_parities[i])
      else:
        hedge_factor *= dual_C_coeff(G.M[edge_list[i][1],edge_list[i][0]][1],G.M[edge_list[i][2],edge_list[i][0]][1],edge_parities[i])
    hedge_factors[parity] += hedge_factor
  return hedge_factors

@cached_function  
def FZ_kappa_factor_graph(G,sigma):#(num,sigma,g,r,markings=(),moduli_type=MODULI_ST):
  #global FZ_kappa_factor_dict
  #key = (num,sigma,g,r,markings)
  #if not FZ_kappa_factor_dict.has_key(key):
  #  G = get_single_stratum(num,g,r,markings,moduli_type)
    L = []
    nv = G.num_vertices()
    for i in range(1,nv+1):
      L.append((2*G.M[i,0][0]+G.degree(i)-2,G.M[i,0]-G.M[i,0][0]))
    LL = []
    tau = []
    for i in range(nv):
      min = -1
      for j in range(nv):
        if (i == 0 or L[j] > LL[-1] or L[j] == LL[-1] and j > tau[-1]) and (min == -1 or L[j] < L[min]):
          min = j
      tau.append(min)
      LL.append(L[min])
    factor_dict = FZ_kappa_factor(tuple(LL),sigma)#get_FZ_kappa_factor2(tuple(LL),sigma)
    factor_vec = [0 for i in range(1 << nv)]
    for parity_key in factor_dict.keys():
      parity = 0
      for i in range(nv):
        if parity_key[i] == 1:
          parity += 1 << tau[i]
      factor_vec[parity] = factor_dict[parity_key]
  #  FZ_kappa_factor_dict[key] = factor_vec
  #return FZ_kappa_factor_dict[key]
    return factor_vec

@cached_function    
def FZ_kappa_factor(L,sigma):
  nv = len(L)
  mmm = max((0,)+sigma)
  sigma_grouped = [0 for i in range(mmm)]
  for i in sigma:
    sigma_grouped[i-1] += 1
  S_list = []
  for i in sigma_grouped:
    S_list.append(IntegerVectors(i,nv))
  S = itertools.product(*S_list)
  kappa_factors = {}
  for parity in itertools.product(*((0,1) for i in range(nv))):
    kappa_factors[tuple(parity)] = 0
  for assignment in S:
    assigned_sigma = [[] for j in range(nv)]
    for i in range(mmm):
      for j in range(nv):
        for k in range(assignment[i][j]):
          assigned_sigma[j].append(i+1)
    sigma_auts = Integer(1)
    parity = [0 for i in range(nv)]
    kappa_factor = Integer(1)
    for j in range(nv):
      sigma_auts *= aut(assigned_sigma[j])
      parity[j] += sum(assigned_sigma[j])
      parity[j] %= 2
      kappa_factor *= kappa_coeff(assigned_sigma[j],L[j][0],L[j][1])
    kappa_factors[tuple(parity)] += kappa_factor/sigma_auts
  return kappa_factors
  
@cached_function
def FZ_marking_factor(G,marking_vec):
  nv = G.num_vertices()
  ne = G.num_edges()
  num_parities = 2**nv
  PPP_list = []
  for marks in marking_vec:
    PPP_list.append(Permutations(marks[1]))
  PPP = itertools.product(*PPP_list)
  marking_factors = [0 for i in range(num_parities)]
  incident_vertices = []
  for mark_type in marking_vec:
    incident_vertices.append([])
    for k in range(1,ne+1):
      if G.M[0,k] == mark_type[0]:
        for i in range(1,nv+1):
          if G.M[i,k] != 0:
            incident_vertices[-1].append((i-1,G.M[i,k][1]))
            break
  for perms in PPP:
    parity = 0
    marking_factor = 1
    for marks_index in range(len(marking_vec)):
      for count in range(len(incident_vertices[marks_index])):
        marking_factor *= C_coeff(perms[marks_index][count],incident_vertices[marks_index][count][1])
        parity ^= (perms[marks_index][count] % 2) << incident_vertices[marks_index][count][0] #should be xor
    marking_factors[parity] += marking_factor
  return marking_factors

@cached_function(key = lambda i,j,parity : (i,j,parity%2))
def dual_C_coeff(i,j,parity):
  #global dual_C_coeff_dict
  #key = (i,j,parity % 2)
  #if dual_C_coeff_dict.has_key(key):
  #  return dual_C_coeff_dict[key]
  total = 0
  k = parity % 2
  while (floor(k/3) <= i):
    if (k % 3) == 2:
      k += 2
      continue
    total += (-1)**(floor(k/3))*C_coeff(k,i)*C_coeff(-2-k,j)
    k += 2
  #dual_C_coeff_dict[key] = total
  return total
  

        
def degenerate1_list(G):
    """
    Returns the degenerations of a single graph. 
    Pasted and edited from tautrel.
    """
    G_list_new = []
    for i in range(1,G.num_vertices()+1):
        #print "i", i
        vi_genus = G.M[i,0]
        if vi_genus > 0:
            G_copy = StrataGraph(G.M)
            G_copy.add_edge(i,i) #add a loop
            G_copy.M[i,0] -= 1 #reduce the genus
            G_copy.compute_invariant()
            G_list_new.append(G_copy) 
        
        row = list(G.M[i])
        m = row[0] + sum(row)
        if m < 4:
            continue
        row1 = [0 for j in range(len(row))]
        #print "just before while"
        while [2*x for x in row1] <= row:
            if row1[0] + sum(row1) >= 2 and row1[0] + sum(row1) <= m-2:
                row2 = [row[j] - row1[j] for j in range(len(row))]
                G_copy = StrataGraph(G.M)
                G_copy.split_vertex(i,row1,row2)
                G_copy.compute_invariant()
                G_list_new.append(G_copy)
            row1[-1] += 1
            for j in range(1,len(row)):
                if row1[-j] <= row[-j]:
                    break
                row1[-j] = 0
                row1[-j-1] += 1
          #print "finished a while loop", [2*x for x in row1], row
    #print "ready to remove dupes"

    return G_list_new

def dim_form(g,n):
  return 3*g-3+n
      
def decorate1_list(G,r):
    X = StrataGraph.Xvar
    G_deco = []
    G.compute_degree_vec()
      #print "compute_degree_vec ok"
    nr,nc = G.M.nrows(),G.M.ncols()
    two_list = []
    one_list = []
    for i in range(1,nr):
        for j in range(1,nc):
            if G.M[i,j] == 2:
                two_list.append([i,j])
            elif G.M[i,j] == 1:
                one_list.append([i,j])
    a = nr-1
    b = len(two_list)
    c = len(one_list)

    #dims = [[dim_form(G.M[i+1,0][0], G.degree(i+1), mod_type) for i in range(a)] for mod_type in range(moduli_type+1)]
    dims = [dim_form(G.M[i+1,0][0], G.degree(i+1)) for i in range(a)] 

    for vec in IntegerVectors(r,a+b+c):
        bad = False
        test_dims = vec[:a]
        for i in range(b):
            test_dims[two_list[i][0]-1] += vec[a+i]
        for i in range(c):
            test_dims[one_list[i][0]-1] += vec[a+b+i]
        #new_type = which_type
        #for mod_type in range(which_type,moduli_type+1):
        for i in range(a):
            if test_dims[i] > dims[i]:
                bad = True
                break
        #            new_type = mod_type + 1
        #            break
        if bad:
            continue
        #if new_type > moduli_type:
        #    continue
        S_list = []
        for i in range(a):
            #print i, vec, vec[i] 
            S_list.append(Partitions(vec[i]))
        for i in range(a,a+b):
            S_list.append([[vec[i]-j,j] for j in range(Rational((vec[i],2)) + 1)])
        S = itertools.product(*S_list) #removed a * here
        
        for vec2 in S:
            G_copy = StrataGraph(G.M)
            for i in range(a):
                for j in vec2[i]:
                    G_copy.M[i+1,0] += X**j
            for i in range(a,a+b):
                G_copy.M[two_list[i-a][0],two_list[i-a][1]] += vec2[i][0]*X + vec2[i][1]*X**2
            for i in range(c):
                G_copy.M[one_list[i][0],one_list[i][1]] += vec[i+a+b]*X
            G_copy.compute_invariant()
            G_deco.append(G_copy) #[new_type].append(G_copy)
      #print "about to return"
    return G_deco
    
   
def kappa_coeff_key(sigma,kappa_0,F):
  mmm = F.degree()
  target_partition = []
  for i in range(1,mmm+1):
    for j in range(F[i]):
      target_partition.append(i)
  key = (tuple(sigma),kappa_0,tuple(target_partition))
  return key

@cached_function(key = kappa_coeff_key)
def kappa_coeff(sigma,kappa_0,F):
  #maybe not optimal to compute the target_partition twice here...
  #global kappa_coeff_dict
  mmm = F.degree()
  target_partition = []
  for i in range(1,mmm+1):
    for j in range(F[i]):
      target_partition.append(i)
  #key = (tuple(sigma),kappa_0,tuple(target_partition))
  #if kappa_coeff_dict.has_key(key):
  #  return kappa_coeff_dict[key]
  total = 0
  num_ones = sum(1 for i in sigma if i == 1)
  for i in range(0,num_ones+1):
    for injection in Permutations(range(len(target_partition)),len(sigma)-i):
      term = binomial(num_ones,i)*binomial(kappa_0 + len(target_partition) + i-1, i)*factorial(i)
      for j in range(len(sigma)-i):
        term *= C_coeff(sigma[j+i],target_partition[injection[j]])
      for j in range(len(target_partition)):
        if j in injection:
          continue
        term *= C_coeff(0,target_partition[j])
      total += term
  return (-1)**(len(target_partition)+len(sigma))*total * Rational((1,aut(target_partition)))
  #kappa_coeff_dict[key] = (-1)**(len(target_partition)+len(sigma))*total/aut(target_partition)
  #return kappa_coeff_dict[key]

@cached_function  
def reduced_FZ_param_list(G,v,g,d,n):
  X = StrataGraph.Xvar
  params = FZ_param_list(n,tuple(range(1,d+1)))
  graph_params = []
  M = Matrix(StrataGraph.R,2,d+1)
  M[0,0] = -1
  for i in range(1,d+1):
    M[0,i] = i
  for p in params:
    G_copy = StrataGraph(G.M)
    M[1,0] = -g-1
    for j in p[0]:
      M[1,0] += X**j
    for i in range(1,d+1):
      M[1,i] = 1 + p[1][i-1][1][0]*X
    G_p = StrataGraph(M)
    G_copy.replace_vertex_with_graph(v,G_p)
    graph_params.append([p,G_copy])
  params_reduced = []
  graphs_seen = []
  for x in graph_params:
    x[1].compute_invariant()
    good = True
    for GG in graphs_seen:
      if graph_isomorphic(x[1],GG):
        good = False
        break
    if good:
      graphs_seen.append(x[1])
      params_reduced.append(x[0])
  return params_reduced  
  
@cached_function   
def FZ_param_list(n,markings):
    #n= 3*r-self.genus-1 #this is how I've seen it called
    #markings = self.markings
    
    if n < 0:
        return []
    final_list = []
    mmm = max((0,)+markings)
    markings_grouped = [0 for i in range(mmm)]
    for i in markings:
        markings_grouped[i-1] += 1
    markings_best = []
    for i in range(mmm):
        if markings_grouped[i] > 0:
            markings_best.append([i+1,markings_grouped[i]])
    for j in range(n/2 + 1):
        for n_vec in IntegerVectors(n-2*j,1+len(markings_best)):
            S_list = [[list(sigma) for sigma in Partitions(n_vec[0]).list() if sum(1 for l in sigma if (l % 3) == 2) == 0]]
            for i in range(len(markings_best)):
                S_list.append(Partitions(n_vec[i+1]+markings_best[i][1], length=markings_best[i][1]).list())
                S_list[-1] = [[k - 1 for k in sigma] for sigma in S_list[-1] if sum(1 for l in sigma if (l % 3) == 0) == 0]          

            for S in itertools.product(*S_list):
                final_list.append((tuple(S[0]), tuple([(markings_best[k][0], tuple(S[k+1])) for k in range(len(markings_best))])))
    return final_list

  
def C_coeff(m,term):
  n = term - floor(m/3)
  if n < 0:
    return 0
  if (m % 3) == 0:
    return A_list[n]
  else:
    return B_list[n]
    
def tuplify(S_list):
    """
    Drew added this.
    This should change S_list from a list of lists of lists to a list of lists of tuples.
    """
    raise Exception("depreciated!!!")
    return [[tuple(l2)  for l2 in l1] for l1 in S_list]