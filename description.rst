Installation
==============

The module is distributed on PyPI. You just need to run the following command: ::

    $ sage -pip install mgn --user

The ``–-user`` option is optional and allows you to install the module in your user space (and does not require administrator rights). 

If this doesn't work (perhaps because your installation of Sage doesn't support SSL), click on the download link and extract the files. Open a terminal and go to the extracted directory ``mgn`` (which contains the ``setup.py`` file). Run: ::

    $ sage -pip install .

mgn Package
============

The ``mgn`` package contains two modules for computations on the moduli space of curves: ``strataalgebra`` and ``topintersections``.

strataalgebra
--------------

This module computes product in the strata algebra. It also contains some of A. Pixton's code to compute FZ relations. `Read the documentation <https://rawgit.com/uberparagon/mgn/master/strataalgebra/_build/html/index.html>`__.


To see if it is installed, try: ::
    
    sage: from strataalgebra import *
    sage: StrataAlgebra(QQ, 1, (1,2))
    Strata algebra with genus 1 and markings (1, 2) over Rational Field
    
topintersections
----------------- 

This module computes top intersections. `Read the documentation <https://rawgit.com/uberparagon/mgn/master/topintersections/_build/html/index.html>`__.

To see if it is installed, try: ::
    
    sage: from topintersections import *
    sage: space(1, 2)
    [1]  psi1
    [2]  psi2
    [3]  ka1
    [4]  ka2
    [5]  ch1
    [6]  irr
    [7]  Dg0m1_2
    [8]  la1
    Mbar_1_2
    sage: intnum(psi1*psi2)
    Computing the intersection of psi1*psi2 over Mbar_1_2...
    1/24
