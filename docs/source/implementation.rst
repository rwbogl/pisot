Maple and Python Comparison
===========================

Doron Zeilberger's original `Pisot.txt`_ was implemented in the computer
algebra system Maple. Maple has a very robust symbolic computation system, but
its programming language is clunky. I have chosen to implement `Pisot.txt`_ in
the programming language Python to enhance its readability and show that Maple
(and its ilk) are not the only languages that can perform symbolic computation.

Python was designed to be a readable and general purpose language. Applications
specific to mathematicians are not very general purpose, so symbolic
computation is not a native feature of Python. Fortunately the library `SymPy`_
implements symbolic computation in Python and along with it most features found
in Maple. By performing all of the "behind the scenes" calculation with SymPy,
we obtain both the readability of Python and the symbolic power of Maple.

.. _Pisot.txt: http://sites.math.rutgers.edu/~zeilberg/tokhniot/Pisot.txt

.. _SymPy:

As a short example of the differences between Maple and Python, we will compare
two snippets from Zeilberger's Maple procedure ``Pis`` and the corresponding
Python version. The procedure ``Pis`` computes the absolute value of the second
largest root of the characteristic equation for a C-finite recurrence. As a
part of this procedure, Zeilberger discards ``1`` if it is a root. If ``lu``
contains the roots, then the following Maple snippet accomplishes this task::

    if member(1,lu) then
    lu:=convert({op(lu)} minus {1},list):
    if lu=[] then
     RETURN(FAIL):
    fi:
    fi:

Perhaps this is an obvious design pattern to an experienced Maple programmer.
To my eyes, there is a cognitive barrier between what is being done (removing
the root 1) and what is being written (convert lists to sets and doing set
subtraction and then converting back to lists).

If ``roots`` is a dictionary of roots, then in Python the same task can be
accomplished as::

    if 1 in roots:
        del roots[1]

    if not roots:
        return None

The Python version communicates the intent of the code more clearly and without
as many technical details.

Detailed Comparison
-------------------

As an example of the differences between the two, let us compare the full
implementations of Zeilberger's Maple procedure ``Pis``. His procedure is as
follows::

    #Pis(C): Inputs a C-finite sequence and outputs the absolute value of the second-largest root
    #It is a Pisot number if it is less than 1.
    #Fis([[1,1],[1,1]]);
    #Pis([[10,219,4796,105030],[22,-3,18,-11]]);
    Pis:=proc(C) local x,lu,i,aluf,mu:
    lu:=[solve(x^nops(C[2])-add(C[2][i]*x^(nops(C[2])-i),i=1..nops(C[2])))]:

    if nops(lu)<>nops(C[1]) then
     RETURN(FAIL):
    fi:

    if member(1,lu) then
    lu:=convert({op(lu)} minus {1},list):
    if lu=[] then
     RETURN(FAIL):
    fi:
    fi:

    aluf:=1:

    for i from 2 to nops(lu) do
     if abs(evalf(lu[i]))>abs(evalf(lu[aluf])) then
      aluf:=i:
     fi:
    od:

    mu:=evalf([op(1..aluf-1,lu),op(aluf+1..nops(lu),lu)]):

    max(seq(abs(mu[i]),i=1..nops(mu))):


    end:

Now, taking our implementation of :class:`cfinite.CFinite` for granted, the
Python
implementation is as follows::

    def pisot_root(c_seq):
        """
        Compute the absolute value of the second-largest root of the characteristic
        equation for the C-finite sequence.

        :c_seq: :class:`.CFinite` instance.

        :returns: Floating point evaluation of the absolute value of the root, or
                  None.

        """
        roots = c_seq.characteristic_roots()
        n_roots = len(roots.keys())

        if n_roots != c_seq.degree:
            return None

        if 1 in roots:
            del roots[1]

        if not roots:
            return None

        root_norms = [abs(root) for root in roots.keys()]
        root_norms = [sympy.re(sympy.N(norm)) for norm in root_norms]

        max_index = root_norms.index(max(root_norms))
        del root_norms[max_index]

        return max(root_norms)

The procedures are about the same length, accounting for blank lines and
comments.

The first feature is the documentation. Both versions document their inputs and
outputs, but the Python version is written in a standard way that allows for
automatic documentation generation. In fact, the documentation at
:func:`pisot.pisot_root` is automatically generated, as is every other piece of
documentation on this site.

The next part computes the roots of the characteristic polynomial. Maple::

    lu:=[solve(x^nops(C[2])-add(C[2][i]*x^(nops(C[2])-i),i=1..nops(C[2])))]:

Python (taking advantage of :meth:`cfinite.CFinite.characteristic_roots`)::

    roots = c_seq.characteristic_roots()

Using Python's classes, it is clear that the characteristic roots are a
property of the C-finite sequence, and that these roots are what we are
computing.

Next, we look to see if there are any repeated roots. This is true if the
number of distinct roots is less than the degree of the sequence. Maple::

    if nops(lu)<>nops(C[1]) then
     RETURN(FAIL):
    fi:

Python::

    n_roots = len(roots.keys())

    if n_roots != c_seq.degree:
        return None

The Maple version counts the number of coefficients in the C-finite sequence
and calls this the degree. The structure of the Maple C-finite sequences is
``[[coeffs], [initials]]``, so it counts the number of elements in the first
element of a nested list. The Python version, relying on the class
:class:`cfinite.CFinite`, simply asks for the sequence's degree. Again,
Python's classes allow us to embed information into an object, rather than
relying on its actual representation.
