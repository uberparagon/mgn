"""
This file implements the user interface.  This is the file that should be loaded by the user.  It will import everything you need.
"""

try:
    from intersection7 import *
    from TautRing3 import Mgn
    from sage.all import *
except ImportError:
    pass


default_file = "mgn_top_intersect.dat"
import sys
#if len(sys.argv) >= 2 and sys.argv[1] == "-i":        #Print this message if you load from interactive sage.
if False: #don't print this message any more
    print """***************************************************************************
Welcome to the intersection number computer for the moduli space of curves!

Written by Drew Johnson, based on Carel Faber's exposition and Maple code.
 
The main commands are "intnum" to compute the intersection number, and 
"space" to select the space.  Type help(intnum) or help(space) for more 
information.
 
Type save_data("filename") to save the results computed this session and 
load_data("filname") to retrieve previously saved data.  The default filename is 
"{0}" if not specified.

Ctrl-D to quit.
***************************************************************************""".format(default_file)

current_space = None
def space(g,n, namespace = None, print_classes = True):
    """
    INPUT:
     - ``g`` -- The genus.
     - ``n`` -- The number of marked points.
     - ``print_classes`` -- (optional, defaults to True) Set this to False if you don't want to see a list of the classes    
     
    This performs three functions:
        1.  It prints a list of classes and indexes for easy reference.
        2.  It sets the defaut space for ``intnum``.
        3.  It injects the variables names into your namespace so you can use them to make polynomials.
    
    It also returns a Moduli space object, if you want to capture it.
    """
    global current_space
    current_space = Mgn(g,n,True)
    
    
    if print_classes:
        current_space.rij()           
    #exprwdata.ExprWithData.reset()
    
    return current_space
    
def intnum(*args, **keywrds):
    r"""
    A convenience method for the user to access the intersection number code.  There are several accepted syntaxes.
    
    INPUT:
     - ``p`` -- a polynomial in variables that can be interpreted as classes on `\mathcal M_{g,n}`.  You can inject variables into the global namespace using the ``space(genus, n)`` function. (Use ``space(genus, n, globals())`` if you imported this module instead of loading it.)  The intersection will be computed on the space specified by the most recent ``space`` call.
    
    INPUT:
      - ``l`` -- a list of indexes into the classes on `\mathcal M_{g,n}`.  The intersection will be computed on the space specified by the most recent ``space`` call.
    
    INPUT:
     - ``genus``
     - ``n`` -- number of marked points
     - ``p`` -- a polynomial in variables that can be interpreted as classes on `\mathcal M_{g,n}`.  This method call will interpret the variables on the space you specify, NOT on the space specified by the most recent ``space`` call.  You may have to call ``space``, however, to ensure that the variables you want to use are defined.
     
    INPUT:
     - ``genus``
     - ``n`` -- number of marked points
     - ``l`` -- a list of indexes into the classes on `\mathcal M_{g,n}`.  These indexes are displayed for you when you use the ``space(genus, n)`` command, but it is not necessary to call ``space`` before using this syntax.    
    """
    global current_space
    M = keywrds.get("space")
    if M != None:
        if len(args) != 1:
            print "If you specify the space, you need only one argument!"
        if isinstance(args[0], list):
            p = prod((M[i] for i in args[0]))
        else:
            p = change_space(args[0], M)  
    
    elif len(args) == 1:    
        if isinstance(current_space, type(None)):
            print 'Please specify the genus and number of marked points as the first arguments, or set the default space you wish to work over using the "space" function'
            return
        M = current_space
        if isinstance(args[0], list):
            p = prod((M[i] for i in args[0]))            
        else:
            p = args[0]
        
    elif len(args) == 3:
        M = Mgn(args[0], args[1])
        if isinstance(args[2], list):
            p = prod((M[i] for i in args[2]))
        else:
            p = change_space(args[2],M)
    else:
        print "Syntax incorrect, please see docstring (type ``help(intnum)``)"
        return
    
    if keywrds.get("confirm", True):
        print "Computing the intersection of {0} over {1}...".format(repr(p), repr(M))
    try:
        return intersect([M],p, keywrds.get("check_degree", True))
    except BadDegreeException as excpt:
        print excpt
        
        
    
def change_space(p, space):
    result = 0
    for m, coeff in p.monomials():
        result += coeff * prod([cl.change_space(space)**expon for cl, expon in m.decompose_monomial()])
    return result



    

def save_data(filename = default_file, prompt = True):
    global master_table
    import os
    if prompt and os.path.exists(filename):
        confirm = raw_input("Overwrite existing file " + filename + " (y for yes)?  ")
        if confirm.lower() != "y":
            print "Save aborted."
            return
            
    import cPickle
    try:
        f = open(filename, "wb")
        try:
            cPickle.dump(master_table, f, protocol = 2)
        except Exception as ex:
            print ex
            return
        finally:
            f.close()            
    except IOError:
        print "Error opening file. Could not save." 
    else:
        print "Save suceeded."
    
    
        
def load_data(filename = default_file):
    global master_table
    import cPickle
    try:
        f = open(filename, "rb") 
        try:
            master_table.update(cPickle.load(f))
        except cPickle.UnpicklingError:
            print "Problem loading data... perhaps the data is corrupted or not in the right format?"
            return
        finally:
            f.close()
    except IOError:
        print "Could not load file.  Does the file exist?"
    else:
        print "Data loaded."
