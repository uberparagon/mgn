from __future__ import print_function, absolute_import
import itertools
from sage.all import PolynomialRing, var, factorial, Matrix, ZZ, cached_function, copy, Permutations, SetPartitions, prod, SR, var, sage, Integer, operator
  

class StrataGraph(object):
  R = PolynomialRing(ZZ,"X",1,order='lex')
  Xvar = R.gen()
  
  def __init__(self,M=None,genus_list=None):
    #print genus_list
    if M:
      self.M = copy(M)
    elif genus_list:
      self.M = Matrix(StrataGraph.R,len(genus_list)+1,1,[-1]+genus_list)
    else:
      self.M = Matrix(StrataGraph.R,1,1,-1)
      
  def __repr__(self):
    """
    Drew added this.
    """
    name = self.nice_name()
    if name == None:
        return repr(self.nice_matrix())
    else:
        return name

  def num_vertices(self):
    return self.M.nrows() - 1

  def num_edges(self):
    return self.M.ncols() - 1

  def h1(self):
    return self.M.ncols()-self.M.nrows()+1

  def add_vertex(self,g):
    self.M = self.M.stack(Matrix(1,self.M.ncols()))
    self.M[-1,0] = g

  def add_edge(self,i1,i2,marking=0):
    self.M = self.M.augment(Matrix(self.M.nrows(),1))
    self.M[0,-1] = marking
    if i1 > 0:
      self.M[i1,-1] += 1
    if i2 > 0:
      self.M[i2,-1] += 1

  def del_vertex(self,i):
    self.M = self.M[:i].stack(self.M[(i+1):])

  def del_edge(self,i):
    self.M = self.M[0:,:i].augment(self.M[0:,(i+1):])

  def compute_degree_vec(self):
    self.degree_vec = [sum(self.M[i,j][0] for j in range(1,self.M.ncols())) for i in range(1,self.M.nrows())]

  def degree(self,i):
    return self.degree_vec[i-1]

  def split_vertex(self,i,row1,row2):
    self.M = self.M.stack(Matrix(2,self.M.ncols(),row1+row2))
    self.add_edge(self.M.nrows()-2, self.M.nrows()-1)
    self.del_vertex(i)

  def compute_parity_vec_original(self):
    X = StrataGraph.Xvar
    self.parity_vec = [(ZZ(1+self.M[i,0][0]+sum(self.M[i,k][1] for k in range(1,self.M.ncols()))+sum(self.M[i,k][2] for k in range(1,self.M.ncols()))+self.M[i,0].derivative(X).substitute(X=1)) % 2) for i in range(1,self.M.nrows())]
    
  def compute_parity_vec(self):
    X = StrataGraph.Xvar
    self.parity_vec = [(ZZ(1+self.M[i,0][0]+sum(self.M[i,k][1] for k in range(1,self.M.ncols()))+sum(self.M[i,k][2] for k in range(1,self.M.ncols()))+self.last_parity_summand(i)) % 2) for i in range(1,self.M.nrows())]
    
  def last_parity_summand(self,i):
    return sum(( expon[0] * coef for expon, coef in self.M[i,0].dict().items() ))

  def parity(self,i):
    return self.parity_vec[i-1]

  def replace_vertex_with_graph(self,i,G):
    X = StrataGraph.Xvar
    nv = self.num_vertices()
    ne = self.num_edges()
    #i should have degree d, there should be no classes near i, and G should have markings 1,...,d and genus equal to the genus of i
    hedge_list = []
    for k in range(1,self.M.ncols()):
      for j in range(ZZ(self.M[i,k])):
        hedge_list.append(k)
    self.del_vertex(i)
    for j in range(G.num_edges() - len(hedge_list)):
      self.add_edge(0,0)
    for j in range(G.num_vertices()):
      self.add_vertex(G.M[j+1,0])
    col = ne+1
    for k in range(1,G.M.ncols()):
      if G.M[0,k] > 0:
        mark = ZZ(G.M[0,k])
        for j in range(G.num_vertices()):
          if self.M[nv+j,hedge_list[mark-1]] == 0:
            self.M[nv+j,hedge_list[mark-1]] = G.M[j+1,k]
          elif G.M[j+1,k] != 0:
            a = self.M[nv+j,hedge_list[mark-1]][1]
            b = G.M[j+1,k][1]
            self.M[nv+j,hedge_list[mark-1]] = 2 + max(a,b)*X + min(a,b)*X**2
      else:
        for j in range(G.num_vertices()):
          self.M[nv+j,col] = G.M[j+1,k]
        col += 1

  def compute_invariant(self):
    self.compute_parity_vec()
    self.compute_degree_vec()
    nr,nc = self.M.nrows(),self.M.ncols()
    self.invariant = [[self.M[i,0], [], [], [[] for j in range(1,nr)]] for i in range(1,nr)]
    for k in range(1,nc):
      L = [i for i in range(1,nr) if self.M[i,k] != 0]
      if len(L) == 1:
        if self.M[0,k] != 0:
          self.invariant[L[0]-1][2].append([self.M[0,k],self.M[L[0],k]])
        else:
          self.invariant[L[0]-1][1].append(self.M[L[0],k])
      else:
        self.invariant[L[0]-1][3][L[1]-1].append([self.M[L[0],k],self.M[L[1],k]])
        self.invariant[L[1]-1][3][L[0]-1].append([self.M[L[1],k],self.M[L[0],k]])
    for i in range(1,nr):
      self.invariant[i-1][3] = [term for term in self.invariant[i-1][3] if len(term) > 0]
      for term in self.invariant[i-1][3]:
        term.sort()
      self.invariant[i-1][3].sort()
      self.invariant[i-1][2].sort()
      self.invariant[i-1][1].sort()
    vertex_invariants = [[i,self.invariant[i-1]] for i in range(1,nr)]
    self.invariant.sort()
    vertex_invariants.sort(key=lambda x: x[1])
    self.vertex_groupings = []
    for i in range(nr-1):
      if i == 0 or vertex_invariants[i][1] != vertex_invariants[i-1][1]:
        self.vertex_groupings.append([])
      self.vertex_groupings[-1].append(vertex_invariants[i][0])
    
    #Drew added this
    self.invariant = tupleit(self.invariant)
    self.hash = hash(self.invariant)
    
  def __eq__(self,other):
        """
        Drew added this.
        """
        return graph_isomorphic(self, other)
        
  def __hash__(self):
        """
        Drew added this.
        Note that it returns a value stored when "compute_invariant" is called.
        """
        try:
            return self.hash
        except:
            self.compute_invariant()
            return self.hash
            
        
  def codim(self):
    codim = 0
    for v in range(1, self.num_vertices()+1):
        for expon, coef in self.M[v,0].dict().items():
            codim += expon[0]*coef
        for e in range(1, self.num_edges()+1):
            for expon, coef in self.M[v,e].dict().items():
                if expon[0] > 0:
                    codim += coef
    for e in range(1, self.num_edges()+1):
        if self.M[0,e] == 0:
            codim +=1   
    return codim
    
  def codim_undecorated(self):
    codim = 0
    for e in range(1, self.num_edges()+1):
        if self.M[0,e] == 0:
            codim +=1   
    return codim

    print("not implemented!!!!!!")
    return 1

        
  def kappa_on_v(self,v):
    """
    Drew added this.
    
    """
    #print "kappa",self.M[v,0].dict()
    
    for ex, coef in self.M[v,0].dict().items():
        if ex[0] != 0:
            yield ex[0], coef
            
                
  def psi_no_loop_on_v(self,v):
    for edge in range(1,self.num_edges()+1):
        if self.M[v,edge][0] == 1:
            psi_expon = self.M[v,edge][1]
            if psi_expon > 0:
                yield edge, psi_expon
                
  def psi_loop_on_v(self,v):
    for edge in range(1,self.num_edges()+1):
        if self.M[v,edge][0] == 2:
            psi_expon1 = self.M[v,edge][1]
            psi_expon2 = self.M[v,edge][2]
            if psi_expon1 > 0:
                yield edge, psi_expon1, psi_expon2
                    
  def moduli_dim_v(self,v):
    """
    Drew added this.
    """
    return 3*self.M[v,0][0]-3 + sum(( self.M[v,j][0] for j in range(1, self.num_edges()+1 ) ))
   
  def num_loops(self):
    count = 0
    for edge in range(1, self.num_edges()+1):
        for v in range(1, self.num_vertices()+1):
            if self.M[v,edge][0] == 2:
                count+=1
                break
                
    return count
    
  def num_undecorated_loops(self):
    count = 0
    for edge in range(1, self.num_edges()+1):
        for v in range(1, self.num_vertices()+1):
            if self.M[v,edge] == 2:
                count+=1
                break
                
    return count

  def num_full_edges(self):
      count = 0
      for edge in range(1, self.num_edges()+1):
          if sum(( self.M[v,edge][0] for v in range(1, self.num_vertices()+1 ))) == 2:
              count += 1
      return count
    
  def forget_kappas(self):
    M = copy(self.M)
    for v in range(1, self.num_vertices()+1):
        M[v,0] = M[v,0][0]
    return StrataGraph(M)
    
  def forget_decorations(self):
    M = Matrix(StrataGraph.R,[[self.M[r,c][0] for c in range(self.M.ncols())] for r in range(self.M.nrows())] )
    Gnew = StrataGraph(M)
    #Gnew.compute_invariant()
    return Gnew
    
  def is_loop(self,edge):
    for v in range(1, self.num_vertices()+1):
      if self.M[v,edge][0] == 2:
        return True
      if self.M[v,edge][0] == 1:
        return False
        
    raise Exception("Unexpected!")

  ps_name = "ps"
  ps2_name = "ps_"

  def nice_matrix(self):
      Mnice = Matrix(SR, self.M.nrows(), self.M.ncols())
      for edge in range(1, self.num_edges()+1):
          Mnice[0,edge] = self.M[0,edge]
      for v in range(1, self.num_vertices()+1):
          kappas = 1
          for expon, coef in self.M[v,0].dict().items():
              if expon[0]==0:
                  Mnice[v,0] += coef
              else:
                  kappas *= var("ka{0}".format(expon[0]))**coef
          if kappas != 1:
              Mnice[v,0] += kappas

          for edge in range(1, self.num_edges()+1):
              psis = 1
              for expon, coef in self.M[v,edge].dict().items():
                  if expon[0]==0:
                      Mnice[v,edge] += coef
                  elif expon[0]==1:
                      psis *= var(StrataGraph.ps_name)**coef
                  elif expon[0]==2:
                      psis *= var(StrataGraph.ps2_name)**coef
              if psis != 1:
                  Mnice[v,edge] += psis
      return Mnice

  @staticmethod
  def from_nice_matrix(lists):
      Mx = Matrix(StrataGraph.R, len(lists), len(lists[0]))
      Mx[0,0]=-1
      Mx[0,:] = Matrix([lists[0]])
      for v in range(1, len(lists)):
          if lists[v][0] in ZZ:
              Mx[v,0]=lists[v][0]
              continue
          if lists[v][0].operator() == sage.symbolic.operators.add_vararg:
              operands = lists[v][0].operands()
              if len(operands) != 2:
                  raise Exception("Input error!")
              genus = operands[1] #the genus
              kappas = operands[0]
          else:
              kappas = lists[v][0]
              genus = 0


          if kappas.operator() == sage.symbolic.operators.mul_vararg:
              for operand in kappas.operands():
                  Mx[v,0] += StrataGraph._kappaSR_monom_to_X(operand)
              Mx[v,0] += genus
          else:
              Mx[v,0] = StrataGraph._kappaSR_monom_to_X(kappas) + genus

      X = StrataGraph.Xvar
      for v in range(1, len(lists)):
          for edge in range(1, len(lists[0])):
              if lists[v][edge] in ZZ:
                  Mx[v,edge] = lists[v][edge]
              else:
                  for operand in lists[v][edge].operands():
                      if operand in ZZ:
                          Mx[v,edge] += operand
                      elif operand.operator() is None:
                          #it is a ps or ps2
                          Mx[v,edge] += X
                      elif operand.operator() == operator.pow:
                          #it is ps^n or ps2^n
                          Mx[v,edge] += X * operand.operands()[1]
                      elif operand.operator() == sage.symbolic.operators.mul_vararg:
                          #it is a ps^n*ps2^m
                          op1,op2 = operand.operands()
                          if op1.operator() is None:
                              exp1 = 1
                          else:
                              exp1 = op1.operand()[1]
                          if op2.operator() == None:
                              exp2 = 1
                          else:
                              exp2 = op2.operand()[1]

                          if exp1 >= exp2:
                              Mx[v,edge] += exp1*X + exp2*X**2
                          else:
                              Mx[v,edge] += exp1*X**2 + exp2*X
      return StrataGraph(Mx)





  @staticmethod
  def _kappaSR_monom_to_X(expr):
      X = StrataGraph.Xvar
      if expr in ZZ:
          return expr
      elif expr.operator() == None:
          ka_subscript = Integer(str(expr)[2:])
          return  X ** ka_subscript
      elif expr.operator() == operator.pow:
          ops = expr.operands()
          expon = ops[1]
          ka = ops[0]
          ka_subscript = Integer(str(ka)[2:])
          return  expon * X ** ka_subscript

  def nice_name(self):
      """
      :return:
      """
      if self.num_vertices() == 1:
          num_loops = self.num_loops()
          if num_loops >1:
              return None
          elif num_loops == 1:
              if self.codim() == 1:
                  return "D_irr"
              else:
                  return None
          #else, there are no loops
          var_strs = []
          for expon, coef in self.M[1,0].dict().items():
              if expon[0] == 0 or coef == 0:
                  continue
              if coef == 1:
                  var_strs.append("ka{0}".format(expon[0]))
              else:
                  var_strs.append("ka{0}^{1}".format(expon[0],coef))
          for he in range(1,self.num_edges()+1): #should only be half edges now
              if self.M[1,he][1] > 1:
                  var_strs.append("ps{0}^{1}".format(self.M[0,he],self.M[1,he][1]))
              elif self.M[1,he][1] ==  1:
                  var_strs.append("ps{0}".format(self.M[0,he]))
          if len(var_strs) > 0:
              return "*".join(var_strs)
          else:
              return "one"
      if self.num_vertices() == 2 and self.num_full_edges() == 1 and self.codim() == 1:
          #it is a boundary divisor
          v1_marks = [self.M[0,j] for j in range(1, self.num_edges()+1) if self.M[1,j] == 1 and self.M[0,j] != 0]
          v1_marks.sort()
          v2_marks = [self.M[0,j] for j in range(1, self.num_edges()+1) if self.M[2,j] == 1 and self.M[0,j] != 0]
          v2_marks.sort()
          if v1_marks < v2_marks:
              g = self.M[1,0]
          elif v1_marks == v2_marks:
              if self.M[1,0] <= self.M[2,0]:
                  g = self.M[1,0]
              else:
                  g = self.M[2,0]
          else:
              g = self.M[2,0]
              #temp = v1_marks
              v1_marks = v2_marks
              #v2_marks = temp
          if len(v1_marks) == 0:
              return "Dg{0}".format(g)
          else:
              return "Dg{0}m".format(g) + "_".join([str(m) for m in v1_marks])



         
def tupleit(t):
  """
  Drew added this.
  Changes a nested list into a nested tuple.
  Needed so that we can hash the invariant.
  """
  return tuple(map(tupleit, t)) if isinstance(t, list) else t

def graph_isomorphic(G1,G2):
  if "invariant" not in G1.__dict__.keys():
    G1.compute_invariant()
  if "invariant" not in G2.__dict__.keys():
    G2.compute_invariant()
  if G1.invariant != G2.invariant:
    return False
  else:
    return isomorphic(G1.M,G2.M,G1.vertex_groupings,G2.vertex_groupings)

def isomorphic(M1,M2,group1,group2):
  nr,nc = M1.nrows(),M1.ncols()
  PermList = [Permutations(list(range(len(group)))) for group in group1]
  for sigma_data in itertools.product(*PermList):
    sigma = [0 for i in range(nr-1)]
    for i in range(len(group1)):
      for j in range(len(group1[i])):
        sigma[group1[i][j]-1] = group2[i][sigma_data[i][j]]
    good = True
    for i in range(1,nr):
      ii = sigma[i-1]
      for j in range(1,i):
        jj = sigma[j-1]
        L1 = []
        for k in range(1,nc):
          if M1[i,k] != 0 and M1[j,k] != 0:
            L1.append([M1[i,k],M1[j,k]])
        L1.sort()
        L2 = []
        for k in range(1,nc):
          if M2[ii,k] != 0 and M2[jj,k] != 0:
            L2.append([M2[ii,k],M2[jj,k]])
        L2.sort()
        if L1 != L2:
          good = False
          break
      if good == False:
        break     
    if good:
      return True
  return False

@cached_function
def graph_count_automorphisms(G,vertex_orbits=False):
  return count_automorphisms(G.M,G.vertex_groupings,vertex_orbits)

def count_automorphisms(M,grouping,vertex_orbits=False):
  nr,nc = M.nrows(),M.ncols()
  count = 0
  PermList = [Permutations(list(range(len(group)))) for group in grouping]
  if vertex_orbits:
    isom_list = []
  for sigma_data in itertools.product(*PermList):
    sigma = [0 for i in range(nr-1)]
    for i in range(len(grouping)):
      for j in range(len(grouping[i])):
        sigma[grouping[i][j]-1] = grouping[i][sigma_data[i][j]]
    good = True
    for i in range(1,nr):
      ii = sigma[i-1]
      for j in range(1,i):
        jj = sigma[j-1]
        L1 = []
        for k in range(1,nc):
          if M[i,k] != 0 and M[j,k] != 0:
            L1.append([M[i,k],M[j,k]])
        L1.sort()
        L2 = []
        for k in range(1,nc):
          if M[ii,k] != 0 and M[jj,k] != 0:
            L2.append([M[ii,k],M[jj,k]])
        L2.sort()
        if L1 != L2:
          good = False
          break
      if good == False:
        break     
    if good:
      count += 1
      if vertex_orbits:
        isom_list.append(sigma)

  if vertex_orbits:
    orbit_list = []
    vertices_used = []
    while len(vertices_used) < nr-1:
      i = [ii for ii in range(1,nr) if ii not in vertices_used][0]
      orbit = []
      for sigma in isom_list:
        if sigma[i-1] not in orbit:
          orbit.append(sigma[i-1])
          vertices_used.append(sigma[i-1])
      orbit.sort()
      orbit_list.append(orbit)
    return orbit_list

  for i in range(1,nr):
    for k in range(1,nc):
      if M[i,k][0] == 2 and M[i,k][1] == M[i,k][2]:
        count *= 2
    L = []
    for k in range(1,nc):
      if M[i,k] != 0:
        if sum(1 for j in range(1,nr) if M[j,k] != 0) == 1:
          L.append([M[0,k],M[i,k]])
    count *= aut(L)

    for j in range(1,i):
      L = []
      for k in range(1,nc):
        if M[i,k] != 0 and M[j,k] != 0:
          L.append([M[i,k],M[j,k]])
      count *= aut(L)
  return count

def aut(L):
  if len(L) == 0:
    return 1
  L.sort()
  total = 1
  n = 1
  last = L[0]
  for l in L[1:]:
    if l == last:
      n += 1
    else:
      n = 1
    total *= n
    last = l
  return total


def PsiKappaVars(A, kappa_only = False):  
    kappa_names = ["kappa{0}_{1}".format(a,v) for v in range(1,A.num_vertices()+1) for a in range(1,A.moduli_dim_v(v)+1)]
    
    if kappa_only:
        PsiKappaRing = PolynomialRing(ZZ,  kappa_names, len(kappa_names) )
    else:
        psi_names = ["psi{0}_{1}".format(e,v) for v in range(1,A.num_vertices()+1) for e in range(1,A.num_edges()+1) if A.M[v,e][0] > 0]
    
        psi_names += ["psi{0}_{1}_2".format(e,v) for v in range(1,A.num_vertices()+1) for e in range(1,A.num_edges()+1) if A.M[v,e][0] == 2] 
    
        PsiKappaRing = PolynomialRing(ZZ,  kappa_names +psi_names)
        
    kappa = {(a,v) : PsiKappaRing.gens_dict()["kappa{0}_{1}".format(a,v)] for v in range(1,A.num_vertices()+1) for a in range(1,A.moduli_dim_v(v)+1) } 
    
    #print PsiKappaRing
    if kappa_only:
        return PsiKappaRing, kappa
    else:
        psi = {(e,v) : PsiKappaRing.gens_dict()["psi{0}_{1}".format(e,v)] for v in range(1,A.num_vertices()+1) for e in range(1,A.num_edges()+1) if A.M[v,e][0] > 0 }
    
        psi2 = {(e,v) : PsiKappaRing.gens_dict()["psi{0}_{1}_2".format(e,v)] for v in range(1,A.num_vertices()+1) for e in range(1,A.num_edges()+1) if A.M[v,e][0] == 2 }
        
        return PsiKappaRing, kappa, psi, psi2
        

def X_to_kappa(f,kappa):
    """
    f -- A polynomial in X. The constant term is ignored.
    kappa -- A function so that kappa(a) is kappa_a.
    
    Returns::
        A polynomial in kappas, converting ...
    """
    psi_list = []
    for expon,coef in f.dict().items():
        #print expon[0], coef
        if expon[0] !=0:
            psi_list += [expon[0]]*coef
        
    #print psi_list
    result = 0 
    for p in SetPartitions(list(range(len(psi_list)))):   
        #print p  
        #print [ kappa( sum((psi_list[i] for i in s)) ) for s in p]
        result +=   prod((factorial(len(s) - 1) for s in p)) * prod(( kappa( sum((psi_list[i] for i in s)) ) for s in p ))
    return result
