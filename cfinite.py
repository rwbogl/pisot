#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Translation of (parts of) Doron Zeilberger's CFinite.txt from Maple to Python.

This file contains classes and functions to do two things:

1. Work with linear recurrence relations with constant coefficients, called
C-finite sequences or recurrences.

2. Guess possible C-finite recurrences a given list of terms might satisfy.

"""

import sympy
import itertools
from sympy.series.sequences import SeqBase

class CFinite(SeqBase):

    """

    This class provides procedures for working with linear recurrence relations
    with constant coefficients, called C-finite sequences.

    We inhereit from sympy's SeqBase class, though our use is not currently
    (2017-09-20) optimized for recursion. (This goes for pisot.py as well.)
    SeqBase uses symoy's @cacheit decorator, which we should try to take
    advantage of for recurrences.  (For example, after computing
    CFinite.coeff(100), CFinite.coeff(99) is not cached.)

    """

    def __init__(self, initial, coeffs):
        """Create the CFinite sequence.

        :initial: List of initial terms of length equal to the degree.
        :coeffs: List of coefficients for the sequence, whose length defines
                 the degree. The coefficients should be ordered so that
                 coeffs[0] is the "left-most" coefficient, i.e.,

                    a_n = c_0 a_{n - 1} + ... + c_{d - 1} a_{n - d}.

        """
        if len(initial) != len(coeffs):
            raise ValueError("Number of initial terms must equal degree")

        self.initial = initial
        self.coeffs = list(coeffs)
        self.degree = len(coeffs)

    @property
    def start(self):
        return 0

    @property
    def stop(self):
        return sympy.oo

    def _eval_coeff(self, index):
        return self.get_terms(index + 1)[-1]

    def get_terms(self, k):
        """Compute the first k terms as a list."""
        return list(itertools.islice(self.gen(), k))

    def gen(self):
        """Yield a generator of terms.

        To compute terms, we will need to keep track of the previous `degree`
        terms. We need to access all of them, but also constantly delete the
        first term and add a new one. For simplicity, we currently use a list
        for this, despite the O(n) deletion. It is possible to use cyclic lists
        for O(1) runtime in both, should we wish later.

        """
        for term in self.initial:
            yield term

        prevs = list(self.initial)

        while True:
            # prev[k] is the term a_{n - d + k} = a_{n - (d - k)}, so we need to
            # multiply it by c_{d - k - 1}.
            next = sum(prevs[k] * self.coeffs[self.degree - k - 1]
                        for k in range(self.degree))
            yield next

            del prevs[0]
            prevs.append(next)

    def characteristic_poly(self, var=sympy.symbols("x")):
        """
        Create the characteristic polynomial of the recurrence relation in
        `var`.

        :var: Symbol to use as the polynomial's variable.
        :returns: Sympy expression.

        """
        return (var**self.degree -
                    sum(var ** (self.degree - k - 1) * self.coeffs[k]
                            for k in range(self.degree)))

    def characteristic_roots(self):
        """Compute the roots of the characteristic equation.

        :returns: List of roots returned by :mod:`sympy`.

        """
        return sympy.roots(self.characteristic_poly())

def find_cfinite_recurrence(seq, n_terms):
    """
    Try to guess a C-finite recurrence of the given degree that the first
    n_terms terms might satisfy, using sympy's find_linear_recurrence()
    function.

    :n_terms: Number of terms to check.
    :seq: Sympy sequence.
    :returns: A CFinite instance or None.

    """
    coeffs = seq.find_linear_recurrence(n_terms)

    if coeffs:
        initial = list(itertools.islice(seq.gen(), len(coeffs)))
        return CFinite(initial, coeffs)

    return None

def guess_cfinite_degree(terms, degree):
    """
    Try to guess a C-finite recurrence of the given degree that the terms might
    satisfy.

    If the terms do satisfy a linear recurrence relation, then they satisfy a
    certain linear system, where the coefficients of the recurrence are the
    unknowns, and the terms form the coefficient matrix. To try and guess what
    the coefficients should be, we form the system and try to solve it. If it
    works out, then we have a guess. If it doesn't work out, then no possible
    C-finite recurrence can occur of the given degree.

    :terms: Terms of the sequence.
    :degree: Degree of the recurrence to check.
    :returns: A CFinite instance, or None.

    """
    # The given recurrence will require `degree` unknowns to be determined, so
    # we will create `degree` equations to be safe. The first equation requires
    # `degree + 1` terms, the and each one after that adds a new term. Thus, we
    # need `degree + degree = 2 * degree` terms.
    terms_needed = 2 * degree

    if len(terms) < terms_needed:
        err = "The list must be of size >= {} for degree {}".format(terms_needed, degree)
        raise ValueError(err)

    # The equations look like this:
    #       a_{k - d} * c_{d} + ... + a_{k - 1} * c_{1} = a_k.

    rows = []
    for eqn_num in range(degree):
        coeffs = [terms[eqn_num + k] for k in range(degree + 1)]
        rows.append(coeffs)

    augmented_matrix = sympy.Matrix(rows)

    # Source: https://stackoverflow.com/a/9493306/6670539
    symbols = sympy.symbols("c1:{}".format(degree + 1))

    soln = sympy.solve_linear_system(augmented_matrix, *symbols)

    if not soln:
        # No possible linear recurrence exists for this degree (or else we
        # would have found it).
        return None

    if len(soln) < degree:
        print("In guessing a recurrence of degree", degree, "for")
        print("\t" + str(terms) + ",")
        print("we have free variables in the linear system.")
        print("The degree is probably too high!")
        print("\nThe discovered solution is")
        print("\t" + str(soln))
        return None

    ordered_soln = [soln[s] for s in symbols]

    # We have a _possible_ recurrence, so let's check that it holds for the
    # remainder of the terms.

    check_nums = len(terms) - terms_needed
    for eqn_num in range(degree, degree + check_nums):
        lhs = sum(terms[eqn_num + k] * ordered_soln[k] for k in range(degree))
        rhs = terms[eqn_num + degree]
        check = sympy.simplify(lhs - rhs)

        if check != 0:
            print("In guessing a recurrence of degree", degree, "for")
            print("\t" + str(terms) + ",")
            print("we discovered a potential recurrence, given by")
            print("\t" + str(soln) + ".")
            print("Unfortunately, this breaks down for a[{}]".format(eqn_num + degree))
            return None

    return CFinite(terms[:degree], ordered_soln)

def guess_cfinite_recurrence(terms, max_degree = None):
    """
    Guess a C-finite recurrence that the terms might satisfy, up to a maximum
    degree.

    Sympy has something that does this more efficiently, but I wanted to write
    my own that gives prettier error messages (see guess_cfinite_degree()).

    :terms: Terms of the sequence.
    :max_degree: Maximum degree to check.
    :returns: A CFinite instance if a guess is found, otherwise None.

    """
    # We need a certain number of terms to even set up a system to guess at a
    # linear recurrence, so make sure we don't try too high of a degree.
    if not max_degree:
        max_degree = len(terms) // 2
    for degree in range(1, max_degree + 1):
        guess = guess_cfinite_degree(terms, degree)

        if guess:
            return guess

    return None
