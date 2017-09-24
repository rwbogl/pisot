.. pisot documentation master file, created by
   sphinx-quickstart on Wed Sep 20 13:50:44 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pisot Documentation
==============================================================================

.. warning::
   These modules are under development. They may be error-prone, give incorrect
   answers, and be generally terrible.

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

Use
---

To use the package, ensure that `SymPy <http://www.sympy.org/en/index.html>`_
is installed, then clone the `GitHub project
<https://github.com/rwbogl/pisot>`_ (via ``git clone``, or downloading a ZIP
file from the page, etc.). On the terminal, navigate to the directory that
``pisot.py`` and ``cfinite.py`` are in, then open a Python terminal. After
this, import the ``pisot`` module with :code:`import pisot` or some variant.
Then every function defined in ``pisot.py`` will be available as
``pisot.func_name``.

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

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
