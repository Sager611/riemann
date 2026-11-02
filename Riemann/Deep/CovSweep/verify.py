#!/usr/bin/env python3
"""Exhaustive covering-radius check for the q-ary codes in codes/*.txt (q = 3 or 4).

File format: '#' comment lines; a header comment gives 'n=<len>', 'R=<radius>' and 'q=<alphabet size>' (q defaults to 4).
Every other non-empty line is a codeword: n symbols from {0,...,q-1}.
Checks: length and alphabet, codewords distinct, and the exact covering radius (every one of the q^n words lies within
Hamming distance R of a codeword).
Method: the set of words within distance k of the code is a bitset over all q^n words (Python int, bit index
sum_j x_j q^j). One step of the Hamming graph changes coordinate j from symbol s to t, i.e. moves the bits whose j-th
digit is s by (t-s)*q^j. Iterate until all q^n bits are set; the number of steps is the covering radius.
Pure Python, no dependencies. Usage: python3 verify.py [files...]
"""
import glob, os, re, sys, time


def digit_masks(q, n):
    """M[j][s] = bitset of the words whose j-th digit equals s."""
    N = q ** n
    full = (1 << N) - 1
    masks = []
    for j in range(n):
        p = q ** j
        row = []
        for s in range(q):
            m, width = ((1 << p) - 1) << (s * p), q * p      # one period of q*p bits
            while width < N:                               # replicate by doubling
                m |= m << width
                width *= 2
            row.append(m & full)
        masks.append(row)
    return masks


def covering_radius(words, q, n):
    N = q ** n
    full = (1 << N) - 1
    M = digit_masks(q, n)
    S = 0
    for w in words:
        S |= 1 << sum(int(c) * q ** j for j, c in enumerate(w))
    r = 0
    while S != full:
        T = S
        for j in range(n):
            p = q ** j
            for s in range(q):
                part = S & M[j][s]
                if not part:
                    continue
                for t in range(q):
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
    n = R = q = None
    words = []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            for key in ('n', 'R', 'q'):
                m = re.search(r'\b%s=(\d+)' % key, s)
                if m and {'n': n, 'R': R, 'q': q}[key] is None:
                    if key == 'n': n = int(m.group(1))
                    elif key == 'R': R = int(m.group(1))
                    else: q = int(m.group(1))
            continue
        words.append(s)
    q = q or 4
    assert n is not None and R is not None, 'header must give n= and R='
    alpha = set(str(i) for i in range(q))
    for w in words:
        assert len(w) == n and set(w) <= alpha, 'bad codeword %r' % w
    assert len(set(words)) == len(words), 'repeated codeword'
    t = time.time()
    cr = covering_radius(words, q, n)
    ok = cr <= R
    print('%-40s q=%d n=%-2d R=%d |C|=%-6d covering radius=%d  %s  (%.1fs)' % (
        os.path.basename(path), q, n, R, len(words), cr, 'OK' if ok else 'FAIL', time.time() - t))
    return ok


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(here, 'codes', '*.txt')))
    results = [check(f) for f in files]
    print('ALL OK' if results and all(results) else 'SOME FAILED')
    sys.exit(0 if results and all(results) else 1)

