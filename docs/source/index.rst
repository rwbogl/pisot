.. pisot documentation master file, created by
   sphinx-quickstart on Wed Sep 20 13:50:44 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pisot Documentation
==============================================================================

.. toctree::
    :hidden:

    implementation
    modules


Pisot sequences are a family of recursive sequences defined by

.. math::
    \begin{align}
        a_0 &= x \\
        a_1 &= y \\
        a_n &= \left\lfloor \frac{a_{n - 1}^2}{a_{n - 2}} + r \right\rfloor,
    \end{align}

where :math:`0 < x < y` are integers and :math:`r` is some constant. These
modules define procedures to work with Pisot sequences and determine whether or
not they satisfy linear recurrence relations with constant coefficients. This
work is a translation of Doron Zeilberger's Maple package `Pisot.txt
<http://sites.math.rutgers.edu/~zeilberg/tokhniot/Pisot.txt>`_, based on `this
article <https://arxiv.org/abs/1609.05570>`_ by Zeilberger and Neil Sloane. It
is my hope that this Python implementation provides a more accessible example
of the uses of symbolic computation in languages more general than Maple.

To use the package, ensure that `SymPy <http://www.sympy.org/en/index.html>`_
is installed, then clone the `GitHub project
<https://github.com/rwbogl/pisot>`_ (via ``git clone
https://github.com/rwbogl/pisot.git``, or downloading a ZIP file from the page,
etc.). On the terminal, navigate to the directory that ``pisot.py`` and
``cfinite.py`` are in, then open a Python terminal. After this, import the
``pisot`` module with :code:`import pisot` or some variant.  Then every
function defined in ``pisot.py`` will be available as ``pisot.func_name``.

Examples
--------

The highlight of the package is :func:`~pisot.pisot_to_cfinite`. Given a
:class:`~pisot.Pisot` instance, it tries to determine whether or not it
satisfies a linear recurrence relation with constant coefficients.

Example::

    In [1]: import pisot

    In [2]: from sympy import Rational

    In [3]: p = pisot.Pisot(5, 17, Rational(1, 2))

    In [4]: guess_terms = 10

    In [5]: check_terms = 1000

    In [6]: c = pisot.pisot_to_cfinite(p, guess_terms, check_terms)

    In [7]: c
    Out[7]: CFinite([5, 17], [4, -2])

    In [8]: c.get_terms(10)
    Out[8]: [5, 17, 58, 198, 676, 2308, 7880, 26904, 91856, 313616]

    In [9]: p.get_terms(10)
    Out[9]: [5, 17, 58, 198, 676, 2308, 7880, 26904, 91856, 313616]

Or, using the ``verbose`` flag::

    In [10]: pisot.pisot_to_cfinite(p, guess_terms, check_terms, verbose=True)
    The Pisot sequence E_{1/2}(5, 17), whose first few terms are
        [5, 17, 58, 198, 676, 2308, 7880, 26904, 91856, 313616],
    appears to satisfy the C-finite recurrence CFinite([5, 17], [4, -2]) whose first few terms are
        [5, 17, 58, 198, 676, 2308, 7880, 26904, 91856, 313616],

    The absolute value of the second-largest root of the C-finite sequence is 0.585786437626905 <= 1.
    Therefore our conjecture holds.
    Out[10]: CFinite([5, 17], [4, -2])

    In [11]: p = pisot.Pisot(8, 16, Rational(1, 2))

    In [12]: pisot.pisot_to_cfinite(p, guess_terms, check_terms, verbose=True)
    The Pisot sequence E_{1/2}(8, 16), whose first few terms are
        [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096],
    appears to satisfy the C-finite recurrence CFinite([8], [2]) whose first few terms are
        [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096],

    This C-finite sequence looks like a geometric sequence, so this is easier to check.
    The conjecture holds if x divides y; the guessed ratio equals y / x; and r is in [0, 1).

    We already know that r satisfies this.
    The conjectured geometric sequence has ratio 2.

    The ratio is an integer and equals y / x, so our conjecture holds.
    Out[12]: CFinite([8], [2])

Classes
-------

We define the classes :class:`pisot.Pisot` and :class:`cfinite.CFinite`. These
inherit from :class:`.SeqBase`, so they are fully-fledged sequences that can be
used in SymPy. The classes encapsulate some important operations and properties
of C-finite and Pisot sequences, namely:

- Computing lists of terms is done with :meth:`.Pisot.get_terms` and
  :meth:`.CFinite.get_terms`.

- For CFinite sequences, the characteristic polynomial is computed with
  :meth:`.CFinite.characteristic_poly`, and its roots with
  :meth:`.CFinite.characteristic_roots`.

- Guessing a linear recurrence with constant coefficients can be done with
  :meth:`.Pisot.find_cfinite_recurrence`.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
