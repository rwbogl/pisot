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
