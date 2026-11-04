#!/usr/bin/env python3
"""Check of Section 6 (Barker corollary) and small-order sanity checks, by exhaustive bit enumeration.
(1) the paper's table rows are Barker sequences;
(2) all Barker sequences of length 2..20 (first entry +1): lengths found must be {2,3,4,5,7,11,13};
(3) periodic autocorrelations of every sign sequence satisfy P(t) = n mod 4 (n <= 16);
(4) every even-length Barker sequence of length > 2 is a circulant Hadamard row (paper's reduction);
(5) circulant Hadamard rows of order 4 and 16 by exhaustion: 8 at n = 4 (all with row sum +-2), none at n = 16.
"""
def acorr(a, n, t):          # aperiodic, a as bitmask with bit=1 meaning -1
    m = (1 << (n - t)) - 1
    return (n - t) - 2 * (((a ^ (a >> t)) & m).bit_count())
def pcorr(a, n, t):
    full = (1 << n) - 1
    rot = ((a >> t) | (a << (n - t))) & full
    return n - 2 * ((a ^ rot).bit_count())
def from_str(s):
    return sum(1 << k for k, ch in enumerate(s) if ch == '-')

table = {2: '++', 3: '++-', 4: '+++-', 5: '+++-+', 7: '+++--+-', 11: '+++---+--+-', 13: '+++++--++-+-+'}
ok1 = all(all(abs(acorr(from_str(s), n, t)) <= 1 for t in range(1, n)) for n, s in table.items())
print(f"(1) paper's Barker table rows are Barker: {ok1}")

found = {}
ok4 = True
for n in range(2, 21):
    cnt = 0
    for a in range(0, 1 << (n - 1)):          # bit 0 = 0  -> first entry +1
        if all(abs(acorr(a, n, t)) <= 1 for t in range(1, n)):
            cnt += 1
            if n % 2 == 0 and n > 2:
                ok4 &= all(pcorr(a, n, t) == 0 for t in range(1, n))
    if cnt: found[n] = cnt
print(f"(2) Barker sequences (first entry +1) by length 2..20: {found}")
ok3 = all(pcorr(a, n, t) % 4 == n % 4 for n in range(1, 17) for a in range(1 << n) for t in range(1, n))
print(f"(3) P(t) = n mod 4 for all sign sequences, n <= 16: {ok3}")
print(f"(4) even-length Barker (n > 2) rows are circulant Hadamard rows: {ok4}")
res = {}
for n in (4, 16):
    res[n] = [a for a in range(1 << n) if all(pcorr(a, n, t) == 0 for t in range(1, n))]
print(f"(5) circulant Hadamard rows: n=4: {len(res[4])}, n=16: {len(res[16])}")
good = ok1 and set(found) == {2, 3, 4, 5, 7, 11, 13} and ok3 and ok4 and len(res[4]) == 8 and len(res[16]) == 0
print("SUMMARY barker/small orders:", "consistent" if good else "INCONSISTENT")

