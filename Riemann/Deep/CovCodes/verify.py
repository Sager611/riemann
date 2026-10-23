#!/usr/bin/env python3
"""Exhaustive covering-radius check for the quaternary codes in codes/*.txt.

File format: '#' comment lines; one header comment must contain 'n=<len>' and 'R=<radius>'.
Every other non-empty line is a codeword: n symbols from {0,1,2,3}.
Checks: correct length/alphabet, codewords distinct, and the exact covering radius.
Method: the set of words within distance k of the code is a bitset over all 4^n words (Python int, bit index
sum_j x_j 4^j). One step of the Hamming graph changes one coordinate j from symbol s to t, i.e. moves the bits whose
j-th digit is s by (t-s)*4^j. Iterate until all 4^n bits are set; the number of steps is the covering radius.
Pure Python, no dependencies. Usage: python3 verify.py [files...]
"""
import glob, os, re, sys, time


def digit_masks(n):
    """M[j][s] = bitset of words whose j-th digit equals s."""
    N = 4 ** n
    masks = []
    for j in range(n):
        p = 4 ** j
        block = ((1 << p) - 1)                     # p consecutive ones
        row = []
        for s in range(4):
            # pattern: in each period of 4p bits, bits [s*p, (s+1)*p) are set
            period = block << (s * p)
            m, width = period, 4 * p
            while width < N:                       # replicate by doubling
                m |= m << width
                width *= 2
            row.append(m & ((1 << N) - 1))
        masks.append(row)
    return masks


def covering_radius(words, n):
    N = 4 ** n
    full = (1 << N) - 1
    M = digit_masks(n)
    S = 0
    for w in words:
        S |= 1 << sum(int(c) * 4 ** j for j, c in enumerate(w))
    r = 0
    while S != full:
        T = S
        for j in range(n):
            p = 4 ** j
            for s in range(4):
                part = S & M[j][s]
                for t in range(4):
                    if t > s:
                        T |= part << ((t - s) * p)
                    elif t < s:
                        T |= part >> ((s - t) * p)
        r += 1
        if T == S:
            raise ValueError('closure stalled')
        S = T
    return r


def check(path):
    n = R = None
    words = []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            m = re.search(r'\bn=(\d+)', s)
            if m and n is None:
                n = int(m.group(1))
            m = re.search(r'\bR=(\d+)', s)
            if m and R is None:
                R = int(m.group(1))
            continue
        words.append(s)
    assert n is not None and R is not None, 'header must give n= and R='
    for w in words:
        assert len(w) == n and set(w) <= set('0123'), 'bad codeword %r' % w
    assert len(set(words)) == len(words), 'repeated codeword'
    t = time.time()
    cr = covering_radius(words, n)
    ok = cr <= R
    print('%-44s n=%-2d R=%d |C|=%-5d covering radius=%d  %s  (%.1fs)' % (
        os.path.basename(path), n, R, len(words), cr, 'OK' if ok else 'FAIL', time.time() - t))
    return ok


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(here, 'codes', '*.txt')))
    results = [check(f) for f in files]
    print('ALL OK' if results and all(results) else 'SOME FAILED')
    sys.exit(0 if results and all(results) else 1)
