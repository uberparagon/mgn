.. topintersections documentation master file, created by
   sphinx-quickstart on Fri Sep  8 12:55:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

topintersections: Computing top intersections on Mbar_g,n
===========================================================

This page demonstrates how to compute top intersections on Deligne-Mumford compactification of the moduli space of curves using Sage code written by Drew Johnson. A pdf description of the algorithms can be found `here`_. If you have any questions or problems, or if you find the code useful, please contact the author by email: ``werd2.718@gmail.com``.

.. _here: https://github.com/uberparagon/mgn/blob/master/topintersections/tex/Mgn.pdf

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
Acknowledgments
=================

I worked on this project as a graduate student at Brigham Young University while being funded by the NSA grant of my adviser, Tyler Jarvis. Dr. Jarvis also answered many questions and pointed me to helpful resources during this project. The mathematical basis of the algorithms implemented here come from the `paper`_ accompanying Carel Faber's Maple code, and also from Stephanie Yang's `write up`_ about her Macaulay 2 code.  The code computes the Witten tau function using the recursion of `Liu and Xu`_.

.. _paper: http://arxiv.org/abs/alg-geom/9706006
.. _write up: http://j-sag.org/Volume2/jsag-1-2010.pdf
.. _Liu and Xu: http://arxiv.org/abs/0710.5322


Installing and starting
=========================

``topintersections`` is now distributed as part of the ``mgn`` package on PyPI. `Click here`_ for installation instructions.

.. _Click here: https://pypi.python.org/pypi/mgn/
   
Basic operation
================

Start a sage session and import the module: ::

    sage: from topintersections import *

You can use the space command to set the space you are working in. The first argument is the genus, and the last is the number of marked points. The space command does three things:

    1. It prints a list of classes and indexes for easy reference.
    2. It sets the defaut space. This will be the space that the code will work in if you don't specify one in the function calls.
    3. It injects the variables names into your namespace so you can use them to make polynomials. 
    
.. link 

::

    sage: space(2,1)
    [1]  psi1
    [2]  ka1
    [3]  ka2
    [4]  ka3
    [5]  ka4
    [6]  ch1
    [7]  ch3
    [8]  irr
    [9]  Dg1m1
    [10]  la1
    [11]  la2
    Mbar_2_1

The class Dg1m1 is the class corresponding to the reducible boundary divisor where one component has genus 1 and the marked point 1. The class irr corresponds to the class of the irreducible boundary divisor. The classes psi, ka, ch, and la represent psi, kappa, chern character, and lambda classes respectively.

Now we are ready to compute some things. You can type in a polynomial in the classes given by the space command: 

.. link 

::

    sage: intnum(irr^3*psi1)
    Computing the intersection of irr^3*psi1 over Mbar_2_1...
    -11/6
    sage: intnum(ka2^2)
    Computing the intersection of ka2^2 over Mbar_2_1...
    53/5760
    sage: intnum(3*irr^3*psi1 + 6*ka2^2)
    Computing the intersection of 3*irr^3*psi1 + 6*ka2^2 over Mbar_2_1...
    -5227/960

If you are just computing a monomial with no coefficient, you can pass in the indexes as a list. Thus, the following command computes the same number as the first example above. 

.. link 

::

    sage: intnum([8,8,8,1])
    Computing the intersection of irr^3*psi1 over Mbar_2_1...
    -11/6

You can also specify the space you wish to work over in the function call. Any classes that are in the namespace will be interpreted as being in the space you specified. For example: 

.. link 

::

    sage: intnum(2,2, psi1^5)
    Computing the intersection of psi1^5 over Mbar_2_2...
    1/1152
    
The following syntax should be very similar to Carel Faber's MgnLb.txt Maple program: 

.. link 

::

    sage: intnum(2,2,[1,1,1,1,1])
    Computing the intersection of psi1^5 over Mbar_2_2...
    1/1152

However, in order to type in a polynomial in classes, the names must have been created by a previous space command. Thus, if you have not called space with marked points at least 2 in this session, the following will give an error: 

.. link 

::

    sage: intnum(2,2, psi1*psi2^4)
    Traceback (most recent call last):
    ...
    NameError: name 'psi2' is not defined
    
Instead, do something like: 

.. link 

::

    sage: space(2,2)
    [1]  psi1
    [2]  psi2
    [3]  ka1
    [4]  ka2
    [5]  ka3
    [6]  ka4
    [7]  ka5
    [8]  ch1
    [9]  ch3
    [10]  irr
    [11]  Dg0m1_2
    [12]  Dg1m1
    [13]  Dg1m1_2
    [14]  la1
    [15]  la2
    Mbar_2_2
    sage: intnum(psi1*psi2^4)
    Computing the intersection of psi1*psi2^4 over Mbar_2_2...
    1/384

If the degree is not correct, you will know. (The code only computes top intersections.) 

.. link 

::

    sage: intnum(2,2, psi1^2)
    Computing the intersection of psi1^2 over Mbar_2_2...
    The monomial psi1^2 has degree 2, while the space Mbar_2_2 has dimension 5.
    
Some more examples
===================

Here are some of the intersection numbers from Faber's paper: 

.. link 

::

    sage: intnum(4,0, irr^9)
    Computing the intersection of irr^9 over Mbar_4_0...
    -251987683/4320
    sage: intnum(4,0, la1^9)
    Computing the intersection of la1^9 over Mbar_4_0...
    1/113400

Options
========

You can suppress the helpful message that tells you what you are computing using the ``confirm`` keyword argument. 

.. link 

::

    sage: intnum(2,2,psi1^5, confirm = False)
    1/1152

You can have the program return zero if the degree is wrong instead of raising an exception by using the ``check_degree`` keyword argument. 

.. link 

::

    sage: intnum(2,2, psi1^2, check_degree = False)
    Computing the intersection of psi1^2 over Mbar_2_2...

Saving and loading
====================
So far we have computed numbers from scratch. The program automatically saves any answers that it has computed in this session, including numbers computed in recursion steps. For example, if you computed the example ``la1^9`` above, you probably noticed that it took a few seconds. If we compute it again, it will be really fast because the program just looks it up in a dictionary. 

.. link 

::

    sage: timeit("intnum(4,0, la1^9)", number =1, repeat = 1) #random
    1 loops, best of 1: 7.06 ms per loop 

We can save this data to a file to avoid computing it over again in our next session. 

.. link 

::

    sage: save_data("testsave.dat", prompt = False)
    Save suceeded.

If you don't specify a file name, the data is saved to the file ``mgn_top_intersect.dat``.

To load a previously saved data file, use the following command: 

.. link 

::

    sage: load_data("testsave.dat")
    Data loaded.



