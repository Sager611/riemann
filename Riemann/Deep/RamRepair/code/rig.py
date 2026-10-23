# rig.py -- small helpers for outward (directed) conversions between Arb balls and Python numbers.
# Rule used in every script of this directory: a bound is carried as an Arb ball end to end; when a Python double is needed
# (heap keys, loop control, JSON), it is produced by up() / dn(), which round outward: float(x) of an exact Arb number is a faithful
# rounding (one of the two neighbouring doubles), so one further step with math.nextafter lands strictly beyond x.
# Rationals (Fraction) are turned into balls by fr(); ball(lo, hi) is an Arb ball containing the whole interval [lo, hi].
import math
from fractions import Fraction
from flint import arb
def up(a):
    """a double >= every point of the ball a"""
    a = arb(a); u = a.upper()
    if not u.is_finite(): return math.inf
    return math.nextafter(float(u), math.inf)
def dn(a):
    """a double <= every point of the ball a"""
    a = arb(a); l = a.lower()
    if not l.is_finite(): return -math.inf
    return math.nextafter(float(l), -math.inf)
def fr(q):
    """ball containing the rational q (int, Fraction, or a double taken exactly)"""
    q = Fraction(q); return arb(q.numerator)/q.denominator
def ball(lo, hi):
    """ball containing [lo, hi] (lo, hi: doubles, Fractions or Arb balls)"""
    lo = fr(lo) if not isinstance(lo, arb) else lo
    hi = fr(hi) if not isinstance(hi, arb) else hi
    return (lo + hi)/2 + arb(0, ((hi - lo)/2).upper())
def lo_min(*xs):
    """exact lower endpoint of the smallest ball (an Arb number <= every point of every ball)"""
    ls = [arb(x).lower() for x in xs]; m = ls[0]
    for v in ls[1:]:
        if v < m: m = v
    return m
def hi_max(*xs):
    us = [arb(x).upper() for x in xs]; m = us[0]
    for v in us[1:]:
        if v > m: m = v
    return m
def ifloor_up(a):
    """an integer >= floor(x) for every x in the ball a (a >= 0)"""
    return math.floor(up(a))
def s(a, n=30):
    """string of a ball, for JSON; arb(s(a)) contains a"""
    return arb(a).str(n, radius=True)
def upper_of(sv):
    """exact Arb upper endpoint of the ball written as sv"""
    return arb(sv).upper()
def lower_of(sv):
    return arb(sv).lower()
