# pisot

Pisot sequences are a family of recursive sequences involving floors,
introduced by Charles Pisot in 1938. This repository contains Python code to
work with Pisot sequences and determine whether or not they satisfy linear
recurrence relations with constant coefficients. Most of the code is a
translation of Doron Zeilberger's Maple package
[Pisot.txt](http://sites.math.rutgers.edu/~zeilberg/tokhniot/Pisot.txt), based
on [this article](https://arxiv.org/abs/1609.05570) by Zeilberger and Neil
Sloane.

Documentation and details on implementation can be found at the site
http://pisot.readthedocs.io/.

Here is a short example:

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

Or, using the `verbose` flag:

    In [10]: pisot.pisot_to_cfinite(p, guess_terms, check_terms, verbose=True)
    The Pisot sequence E_{1/2}(5, 17), whose first few terms are
        [5, 17, 58, 198, 676, 2308, 7880, 26904, 91856, 313616],
    appears to satisfy the C-finite recurrence CFinite([5, 17], [4, -2]) whose first few terms are
        [5, 17, 58, 198, 676, 2308, 7880, 26904, 91856, 313616],

    The absolute value of the second-largest root of the C-finite sequence is 0.585786437626905 <= 1.
    Therefore our conjecture holds.
