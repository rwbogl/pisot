#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Translation of Doron Zeilberger's Pisot.txt from Maple to Python.

Pisot.txt is a Maple program that implements procedures to assist in the study
of Pisot sequences. In particular, the highlight is its ability to discern (and
prove!) whether or not a given Pisot sequence satisfies a linear recurrence
relation with constant coefficients (up to all degree checked). To do this, it
relies on the Maple program Cfinite.txt, also written by Zeilberger.

Pisot.txt can be difficult to read. Due to its mathematical complexity it is a
fairly large program. This is exacerbated by the fact that Maple, though useful
in its domain, is not an elegant programming language. I hope that this Python
implementation will (eventually) provide a more accessible demonstration of the
applications of symbollic computing.

"""

import itertools
import sympy
from sympy.series.sequences import SeqBase
import cfinite

class Pisot(SeqBase):

    """
    This class defines basic methods for dealing with Pisot sequences.

    The Pisot sequence E_r(x, y) is defined by

        a_0 = x
        a_1 = y
        a_{n + 1} = floor(a_{n - 1}**2 / a_{n - 2} + r).
    """

    def __iter__(self):
        return self.gen()

    def __init__(self, x, y, r):
        """Create the Pisot sequence instance.

        :x: First term in the sequence.
        :y: Second term in the sequence.
        :r: Pisot sequence parameter.

        """
        self.x = x
        self.y = y
        self.r = r

    @property
    def start(self):
        return 0

    @property
    def stop(self):
        return sympy.oo

    def get_terms(self, k):
        """Compute the first k terms as a list."""
        return list(itertools.islice(self.gen(), k))

    def gen(self):
        """Yield a generator of terms."""
        yield self.x
        yield self.y

        back_2 = self.x
        back_1 = self.y

        while True:
            new = sympy.floor(back_1 ** 2 / back_2 + self.r)
            yield new

            back_2, back_1 = back_1, new

    def _eval_coeff(self, index):
        return self.get_terms(index + 1)[-1]

def pisot_root(c_seq):
    """
    Compute the absolute value of the second-largest root of the characteristic
    equation for the C-finite sequence.

    Zeilberger does additional things in this method. If there are repeated
    roots, we return None. If 1 is the only root, we return None.

    :arg c_seq: :class:`.CFinite` instance.

    :returns: Floating point evaluation of the absolute value of the root, or
              None.

    """
    roots = c_seq.characteristic_roots()
    n_roots = sum(1 for root in roots.keys())

    # Repeated roots.
    if n_roots != c_seq.degree:
        return None

    if 1 in roots:
        del roots[1]

    # 1 is the only root.
    if not roots:
        return None

    root_norms = [abs(root) for root in roots.keys()]

    # Sometimes floating point arithmetic doesn't work how we want, so trim off
    # the imaginary parts of the absolute value of a number.
    root_norms = [sympy.re(sympy.N(norm)) for norm in root_norms]

    # Delete the largest root, then return the (next) largest.
    max_index = root_norms.index(max(root_norms))

    del root_norms[max_index]

    return max(root_norms)

def pisot_to_cfinite(pisot, guess_length, check_length):
    """
    Check if the given Pisot sequence satisfies a linear recurrence relation
    with finite coefficients.

    :pisot: :class:`.Pisot` sequence.
    :guess_length: Number of terms to use when guessing the recurrence. This
                   should be somewhat small. If the sequence fails to satisfy a
                   linear recurrence at a large number, then this method will
                   spend a long time trying to look for one.
    :check_length: Number of terms of the sequence to check the conjectured
                    linear recurrence for. This should be a large number.

    :returns: A :class:`.CFinite` instance, or None.

    """
    if pisot.r <= 0 or pisot.r >= 1:
        raise ValueError("r must be strictly between 0 and 1 for this procedure")

    c_seq = cfinite.find_cfinite_recurrence(pisot, guess_length)

    if not c_seq:
        return None

    check_p = pisot.get_terms(check_length)
    check_c = c_seq.get_terms(check_length)

    if check_p != check_c:
        return None

    root_abs = pisot_root(c_seq)

    if root_abs > 1:
        return None

    return c_seq

if __name__ == "__main__":
    sympy.init_printing()
    x = 5
    y = 17
    # Note the use of sympy.Rational here. If we just used (floating point)
    # 1/2, the arithmetic breaks down once things get large.
    r = sympy.Rational(1, 2)
    p = Pisot(x, y, r)
    n_terms = 10
    check_terms = 1000
    c_seq = pisot_to_cfinite(p, n_terms, check_terms)

    # I'm sorry about the triple braces.
    print("The Pisot sequence E_{{{}}}({}, {}),".format(r, x, y), end=" ")
    print("whose first few terms are")
    print("\t" + str(p.get_terms(10)) + ",")
    if c_seq:
        print("satisfies the linear recurrence with constant coefficients")
        print("\t" + str(c_seq.coeffs))
        print("and initial values")
        print("\t" + str(c_seq.initial) + ".")
    else:
        print("does not satisfy any linear recurrence with constant coefficients, up to degree {}.".format(n_terms // 2))
