.. StrataAlgegra documentation master file, created by
   sphinx-quickstart on Thu Aug 31 18:54:11 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

strataalgebra Documentation
=============================
Welcome! The ``strataalgebra`` module is designed to compute products in the strata algebra using `Sage`_. The strata algebra is of interest
because it has a natural surjective morphism to the tautological ring of the moduli space of curves. The kernel is (conjecturally) described by the Faber-Zagier relations. For more detail, see for example `R. Pandharipande's exposition`_.

.. _Sage: sagemath.org
.. _R. Pandharipande's exposition: https://arxiv.org/pdf/1603.05151.pdf

The product structure was implemented by `Drew Johnson`_, based on `S. Yang's note`_ describing the algorithm of Graber and
Pandharipande.
The program also computes the FZ relations using code copied from A. Pixton's `tautrel.sage`_ program.
The code was integrated into this package (with some modifications and optimizations) by Drew Johnson
(who takes responsibility for any bugs introduced!).

.. _Drew Johnson: http://pages.uoregon.edu/drewj/
.. _tautrel.sage: http://math.mit.edu/~apixton/programs/tautrel.sage
.. _S. Yang's note: https://arxiv.org/abs/0808.1974

Installation
=========================================
Installation should be (hopefully) easy. ``strataalgebra`` is distributed as part of the ``mgn`` package on PyPI. `Click here`_ for installation instructions.

.. _Click here: https://pypi.python.org/pypi/mgn/

Once it is installed, you can load it in a sage session by typing: ::

    sage: from strataalgebra import *

How to use
=============

.. autoclass:: strataalgebra.StrataAlgebra
    :members:
    
.. automethod:: strataalgebra.StrataAlgebraElement.integrate
.. automethod:: strataalgebra.StrataAlgebraElement.dict

Testing
=================

.. automodule:: tests

You can also test the examples from this file using Sage's doctest.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

