.. StrataAlgegra documentation master file, created by
   sphinx-quickstart on Thu Aug 31 18:54:11 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

strataalgebra Documentation
=============================
Welcome! The ``strataalgebra`` module is designed to compute products in the strata algebra using sage. The strata algebra is of interest
because it has a natural surjective morphism to the tautological ring of the moduli space of curves. The kernel is described by the Faber-Zagier relations.

The product structure was implemented by `Drew Johnson`_, based on `S. Yang's note`_ describing the algorithm of .....
The program also computes the FZ relations using code copied from A. Pixton's `tautrel.sage`_ program.
The code was integrated into this package (with some modifications and optimizations) by Drew Johnson
(who takes responsibility for any bugs introduced!).

.. _Drew Johnson: http://pages.uoregon.edu/drewj/
.. _tautrel.sage: http://math.mit.edu/~apixton/programs/tautrel.sage
.. _S. Yang's paper: https://arxiv.org/abs/0808.1974

Installation
=========================================
Installation should be (hopefully) easy. ``strataalgebra`` is distributed as part of the ```mgn`` package`_.

.. _``mgn`` package: https://github.com/uberparagon/mgn

Once it is installed, you can load it in a sage session by typing: ::

    sage: from strataalgebra import *

How to use
=============

.. autoclass:: strataalgebra.StrataAlgebra
    :members:

Testing
=================

.. automodule:: tests

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

