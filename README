Installation
==============

The module is distributed on PyPI. You just need to run the following command: ::

    $ sage -pip install surface_dynamics [--user]

The –user option is optional and allows to install the module in your user space (and does not require administrator rights). 

If that doesn't work try updating sage (some earlier versions of sage are packaged with an incompatible version of PyPI.)

If that still fails, try downloading the source, extracting it, opening a terminal in the directory that has the ``setup.py`` file, and run

    $ sage setup.py -install

mgn Package
============

The ``mgn`` package contains two modules for computations on the moduli space of curves.

strataalgebra
==============

This module computes product in the strata algebra. It also contains some of A. Pixton's code to compute FZ relations. 

Read the detailed documentation `here`_.

.. _here: https://rawgit.com/uberparagon/mgn/master/strataalgebra/_build/html/index.html

To see if it is installed, try:

    ::
    
    sage: from strataalgebra import *
    sage: StrataAlgebra(QQ, 1, (1,2))
    StrataAlgebra with genus 1 and markings (1, 2) over RationalField
    
    

topintersections
================= 

To see if it is installed, try:

    ::
    
    sage: from topintersections import *
    ***************************************************************************
    Welcome to the intersection number computer for the moduli space of curves!

    Written by Drew Johnson, based on Carel Faber's exposition and Maple code.
 
    The main commands are "intnum" to compute the intersection number, and 
    "space" to select the space.  Type help(intnum) or help(space) for more 
    information.
 
    Type save_data("filename") to save the results computed this session and 
    load_data("filname") to retrieve previously saved data.  The default filename is 
    "mgn_top_intersect.dat" if not specified.

    Ctrl-D to quit.
    ***************************************************************************