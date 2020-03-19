from __future__ import absolute_import
#load("tautrel.sage")
#load("AllStrata.sage")
from sage.all import Graph, Integer, OrderedSetPartitions, Permutations, flatten
from .stratagraph import StrataGraph
import itertools
   
            
class StrataWithSageGraph(object):
    def __init__(self, strataG):
        self.strataG = strataG
        self.sageG = Graph([list(range(1,strataG.num_vertices()+1)),[]], multiedges=True, loops = True)
        self.edge_label_to_edge = dict()
        self.vertex_to_marks = {v:[] for v in range(1, self.strataG.num_vertices()+1)}
        self._has_marks = False
        
        for e in range(1, strataG.num_edges()+1):
            edge_done = False
            if self.strataG.M[0,e] != 0: #it is a half edge
                for v in range(1, strataG.num_vertices()+1):
                    if self.strataG.M[v,e][0] == 1:
                        self.vertex_to_marks[v].append((e,self.strataG.M[0,e]))
                        edge_done = True
                        self._has_marks = True
                        break
            else: #it is a whole edge           
                   
                vert_touching_e = []
            
                for v in range(1, strataG.num_vertices()+1):
                    if strataG.M[v,e][0] == 2:
                        #add a loop
                        self.sageG.add_edge( (v,v,e) )
                        self.edge_label_to_edge[e] = (v,v,e)
                        edge_done = True
                        break
                    if strataG.M[v,e][0] == 1:
                        vert_touching_e.append(v)
                        if len(vert_touching_e) == 2:
                            break
                            
            if edge_done:
                continue
            if len(vert_touching_e) == 2:
                self.sageG.add_edge( (vert_touching_e[0], vert_touching_e[1],e) )
                self.edge_label_to_edge[e] = (vert_touching_e[0], vert_touching_e[1],e)
            else:
                raise Exception("Unexpected here!")
    
    def has_marks(self):
        return self._has_marks
        
    def marks_on_v(self,v):
        return self.vertex_to_marks[v]
                 
    def edges_incident(self,v):
        for v1,v2,e in self.sageG.edges_incident(v):
            yield e
            
    def edges_labels_between_vertex_sets(self,vs1,vs2):
        #problem here when overlap........
        result = []
        for v1,v2, e in self.sageG.edges():
            if v1 in vs1 and v2 in vs2:
                result.append(e)
            elif v1 in vs2 and v2 in vs1:
                result.append(e)
        return result
        
        #some probably garbage...
        #v1sedges = set()
        #v2sedges = set()
        #for v1 in vs1:
        #    for e in self.edges_incident(v1):
        #        v1sedges.add(e)
        #for v2 in vs2:
        #    for e in self.edges_incident(v2):
        #        v2sedges.add(e)
                
        #return v1sedges.intersection(v2sedges)
        
    def edges(self):
        """
        Returns a list of triples!
        """
        return self.sageG.edges()
        
    def vertices(self):
        return self.sageG.vertices()
        
    def num_vertices(self):
        return self.strataG.num_vertices()
        
    def v_genus(self,v):
        """
        Returns the genus of the vertex v.
        """
        return self.strataG.M[v,0][0]
        
    def subgraph(self,vertices, edge_labels): 
        #print self.edge_label_to_edge
        #print self.sageG.edges()
        return self.sageG.subgraph(vertices, [self.edge_label_to_edge[l] for l in edge_labels])
        
    def edge_is_incident(self,e,v):
        #print "is_inc", e,v
        #print self.strataG
        #print
        return self.strataG.M[v,e][0] > 0
        
    def is_loop(self, e):
        v1,v2, ep = self.edge_label_to_edge[e]
        return v1==v2
                
                

class GStructOnA(object):
    
    def __init__(self,G, A,alpha_inv, beta, gamma_inv):
        self.G = G
        self.A = A
        self.alpha_inv = alpha_inv
        self.beta = beta
        self.gamma_inv = gamma_inv
        
    def __str__(self):
        return """
        ****
        alpha_inv:     {0}
        beta:      {1}
        gamma_inv: {2}
        ****
        """.format(self.alpha_inv, self.beta, self.gamma_inv)
    
    #@cache_method    
    def beta_ray(self,eG,vG):
        """
            vG is a vertex of G.
            eG is an edge of G.
            
            Returns an pair vA,eA
        """
        raise Exception("Depreciated")
        eA = self.beta[eG]
        vA = max((vAi for vAi in self.alpha_inv[vG] if self.A.edge_is_incident(eA,vAi) ))
        return eA, vA
    
    #@cache_method
    def beta_ray2(self,eG,vG):
        raise Exception("Depreciated")
        eA = self.beta[eG]        
        vA = min((vAi for vAi in self.alpha_inv[vG] if self.A.edge_is_incident(eA,vAi) ))
        return eA, vA
        
    def psi_no_loop(self, eG,vG, ex, psi):
        """
        Returns a 
        """
        eA = self.beta[eG]
        vis = [vAi for vAi in self.alpha_inv[vG] if self.A.edge_is_incident(eA,vAi)]
        if len(vis) > 1:
            raise Exception("Unexpected!!!")
        return psi[(eA,vis[0])]**ex
        
    def psi_loop(self, eG, vG, ex1,ex2, psi1,psi2):
        #print self.G.strataG
        #print self.A.strataG
        #print eG
        #print vG
        #print psi1
        #print psi2
        #print "alpha_inv", self.alpha_inv
        #print "beta", self.beta
        eA = self.beta[eG]
        vis = [vAi for vAi in self.alpha_inv[vG] if self.A.edge_is_incident(eA,vAi)]
        #print vis
        #print eA
        #print
        
        if len(vis) == 1:
            return  (psi1[(eA,vis[0])]**ex1 * psi2[(eA,vis[0])]**ex2 + psi1[(eA,vis[0])]**ex2 * psi2[(eA,vis[0])]**ex1)
        if len(vis) == 2:
            return (psi1[(eA,vis[0])]**ex1 * psi1[(eA,vis[1])]**ex2 + psi1[(eA,vis[1])]**ex1 * psi1[(eA,vis[0])]**ex2)
        raise Exception("Unexpected!!!!")
        
    def count_loop_to_loops(self):
        count = 0
        for  vG1, vG2, eG in self.G.edges():
            if vG1 == vG2 and self.A.is_loop(self.beta[eG]):
                count +=1
        return count

# Can't do this until we implement canonical labeling!
#class GStructsOnAMetaclass(type):
#    def __init__(self, *args, **kwargs):
#        super(GStructsOnAMetaclass, self).__init__(*args, **kwargs)
#        self.dict = dict()
#        
#    def __call__(self, *args, **kwargs):
#        G = args[0]
#        A= args[1]
#        
#        value = self.dict.get((G,A),None)
#        if value == None:
#            sp = super(GStructsOnAMetaclass,self).__call__(G,A)
#            self.dict[(G,A)] = sp
#            return sp
#        else:
#            print ".",
#            return value
            
class GStructsOnA(object):
#    __metaclass__ = GStructsOnAMetaclass
    
    def __init__(self,strataG,strataA):
        """
        The inputs should be StrataGraph objects.
        """
        self.G = StrataWithSageGraph(strataG)
        self.A = StrataWithSageGraph(strataA)
        #self._structs = None
        #print "G:"
        #print strataG
        #print self.G.vertex_to_marks
        #print "A:"
        #print strataA
        #print self.A.vertex_to_marks
        
    def __iter__(self):
        for alpha_inv_s in OrderedSetPartitions(self.A.num_vertices(), self.G.num_vertices()):
            
            
            
            beta_marks_list = list(self.get_beta_marks(alpha_inv_s))       
            if len(beta_marks_list) == 0:
                continue
            
            if not self.genus_v_ok(alpha_inv_s):
                continue
                
            #also, we should check the marked points somewhere about here.
            
            #do we need this??
            alpha = dict()    
            for vG,  alpha_inv_vG in zip(self.G.vertices(), alpha_inv_s):
                for vA in alpha_inv_vG:
                    alpha[vA] = vG
                    
            #print "alpha", alpha_inv_s, alpha        
            #print "about to do cp on ",
            #print [list(self.A.edges_labels_between_vertex_sets(alpha_inv_s[v1-1], alpha_inv_s[v2-1]))  for v1,v2,e in self.G.edges() ]
            #print self.G.edges()
                
            
            for beta_s in itertools.product(*[ list(self.A.edges_labels_between_vertex_sets(alpha_inv_s[v1-1],alpha_inv_s[v2-1]))  for v1,v2,e in self.G.edges() ]):
                gamma_inv = {vG:[] for vG in self.G.vertices()}
                #print "beta_s", beta_s
                bad = False
                for vA1, vA2, eA in self.A.edges(): 
                    if not eA in beta_s:
                        if alpha[vA1] != alpha[vA2]:
                            bad = True
                            break
                        else:
                            gamma_inv[alpha[vA1]].append(eA)
                if bad:
                    continue
                        
                        
                        
                            

            
                if not self.genus_ok(alpha_inv_s,gamma_inv):
                    continue
                
                if not self.connected_ok(alpha_inv_s,gamma_inv):
                    continue
                
                #print "about to yield alpha_inv_s" ,alpha_inv_s
                #print "beta_s", beta_s
                alpha_inv = {i+1:alpha_inv_s[i] for i in range(len(alpha_inv_s))}
                beta_no_marks = {eG[2] : eA for eG,eA in zip(self.G.edges(), beta_s)}
                for beta_marks in beta_marks_list:
                    beta = dict(beta_no_marks)
                    beta.update(beta_marks)
                    yield GStructOnA(self.G, self.A, alpha_inv, beta, gamma_inv)
                
    def genus_v_ok(self, alpha_inv_s):
        """
        Checks whether the map of vertices already has too large of genus in the fibers.
        """
        for vG,  alpha_inv_vG in zip(self.G.vertices(), alpha_inv_s):
            if sum(( self.A.v_genus(vA) for vA in alpha_inv_vG)) > self.G.v_genus(vG):
                return False
        return True
        
    def get_beta_marks(self, alpha_inv_s):
        if not self.G.has_marks():
            #print "G has no marks!"
            yield dict()
            return
    
        inv_marks_perms_list = []
        Gmarks = []
        for vG,  alpha_inv_vG in zip(self.G.vertices(), alpha_inv_s):
            inv_marks = []
            for vA in alpha_inv_vG:
                inv_marks += self.A.marks_on_v(vA)
            
            vGmarks = self.G.marks_on_v(vG)    
            Gmarks += vGmarks
            if len(inv_marks) != len(vGmarks):
                return 
            if len(vGmarks) == 0:
                continue
            
            inv_marks_perms = []
            #print "inv_marks", inv_marks
            for p in Permutations(inv_marks):
                bad = False
                for i,j in zip(vGmarks,p):
                    if i[1] != j[1]:
                        bad = True
                        break
                if not bad:
                    inv_marks_perms.append(list(p))
                    #print "inv_marks_perm", inv_marks_perms
            
            if len(inv_marks_perms) > 0:
                inv_marks_perms_list.append(inv_marks_perms)
            else:
                return       
                
        #print "made this",Gmarks, inv_marks_perms_list

                
        for p in itertools.product(*inv_marks_perms_list):
            pf = flatten(p,max_level=1)
            yield {i[0]:j[0] for i,j in zip(Gmarks, pf)}
            
                          
    def genus_ok(self, alpha_inv_s, gamma_inv):
        """
        Checks whether the proposed G-struct on A has the correct genus in the fibers.
        """
        #print alpha_inv_s, gamma_inv
        for v in self.G.vertices():
            inv_genus = sum(( self.A.v_genus(vA) for vA in alpha_inv_s[v-1])) + len(gamma_inv[v]) - len(alpha_inv_s[v-1]) +1
            if inv_genus != self.G.v_genus(v):
                return False
        return True            
        
    def connected_ok(self, alpha_inv_s, gamma_inv):
        """
        Checks whether the fibers of the proposed G-struct on A are connected.
        """
        for v in self.G.vertices():
            invG = self.A.subgraph(alpha_inv_s[v-1], gamma_inv[v])
            if not invG.is_connected():
                return False
        return True

def is_generic(Gs,Hs):
    for v1,v2,e in Gs.A.edges():
        if not (e in Gs.beta.values() or e in Hs.beta.values()):
            return False
    return True
                
def genericGHstructures(G,H,A):
    for Gs in GStructsOnA(G,A):
        for Hs in GStructsOnA(H,A):
            if is_generic(Gs,Hs):
                yield Gs,Hs
                
                

                                
def decorate_with_monomial(PsiKappaRing,A,m):
    """
    A should be a StrataGraph.
    Returns a new StrataGraph decorated according to the monomial m, which is actualy a tuple of exponents for generators of the PsiKappaRing
    """
    #print "entering decorate"
    #print A
    #print m
    X = StrataGraph.Xvar
    Anew = StrataGraph(A.M)
    for gen, e in zip(PsiKappaRing.gens(),m):
        if e == 0:
            continue
        var_str = str(gen)
        if var_str[0] == "k":
            stav = var_str[5:].split("_")
            a= Integer(stav[0])
            v= Integer(stav[1])
            Anew.M[v,0] += e*X**a
        elif var_str[0] == "p":
            stev = var_str[3:].split("_")
            edge = Integer(stev[0])
            v = Integer(stev[1])
            if len(stev) == 2:
                Anew.M[v,edge] += e*X
            if len(stev) == 3:
                if Anew.M[v,edge][1] > 0:
                    if Anew.M[v,edge][1] >= e: #if it already has a psi
                        Anew.M[v,edge] += e*X**2
                    else:
                        Anew.M[v,edge] = Anew.M[v,edge][1]*X**2 + e*X + Anew.M[v,edge][0]
                else: #the psi2 is the only one psi
                    Anew.M[v,edge] += e*X
                        
                
        else:
            raise Exception("You passed the wrong ring!?")
            
    #Anew.compute_invariant()
    #print "returning"
    #print Anew
    return Anew
    
def count_automorphisms_bad(A):
    Ag = StrataWithSageGraph(A)
    return Ag.sageG.automorphism_group().order()
    
def shared_edges(Gs,Hs):
    #print "shared edges for", Gs.beta.values(), Hs.beta.values()
    full_edges = [e for v1,v2,e in Gs.A.edges()]
    edges = set(Gs.beta.values()).intersection(Hs.beta.values()).intersection(full_edges)
    
    for e in edges:
        yield Gs.A.edge_label_to_edge[e]
    
        
#for testing
#Gtest = StrataGraph(Matrix(R,[[-1,0,0,0],[1,2,1,0],[1,0,1,2]]))
#GGtest = StrataWithSageGraph(Gtest)

#for gs in GStructsOnA(Gtest,Gtest):
#    print gs
    
#yG = get_single_stratum(32,4,3,())
#yA = get_single_stratum(855,4,5,())
#yH = get_single_stratum(2185,4,6,())

#results = list(GStructsOnA(yG,yA))

#for r in results:
#    print r
