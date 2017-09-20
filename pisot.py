#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pisot.py - translation of Doron Zeilberger's Pisot.txt from Maple to Python

Pisot.txt is a Maple program that implements procedures to assist in the study
of Pisot sequences. In particular, the highlight is its ability to discern (and
prove!) whether or not a given Pisot sequence satisfies a linear recurrence
relation with constant coefficients (up to all degree checked). To do this, it
relies on the Maple program Cfinite.txt, also written by Zeilberger.

Pisot.txt can be difficult to read. Due to its mathematical complexity, it is a
fairly large program, and this is exacerbated by the fact that Maple, though
useful in its domain, is not an elegant programming language. I hope that this
Python implementation will (eventually) provide a more accessible demonstration
of the applications of symbollic computing.

"""

import itertools
import sympy
from sympy.series.sequences import SeqBase

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

if __name__ == "__main__":
    sympy.init_printing()
    p = Pisot(sympy.E, sympy.pi, 1/2)

    for k, term in enumerate(p.get_terms(6)):
        print("a[{}] = {}".format(k, term))
