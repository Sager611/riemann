#!/usr/bin/env python3
"""Check the 'hardest' production leaves: those containing genuine KKT points (found by kz_locate).
For each: independent exclusion with masses in [0,T] (kz_check) and the rigorous value bound
sup e^{-t}C on the Krawczyk-contracted solution box -> margin to 2/e."""
import sys, json
import kz_check as K
from flint import arb

def value_margin(n, lo, hi):
    N = n - 1
    P = K.domain_contract([[arb(lo[j]), arb(hi[j])] for j in range(N + 1)], N)
    T = arb(K.THI[n])
    Clo, Chi = K.C_range(P, N)
    tcap = (Chi.log() - K.TWO_OVER_E.log()).upper().min(T)
    l = [arb(0)] * N; u = [T] * N
    A, b = K.system_intervals(P, n)
    # without the value cut (to see the solution), certified mass bounds
    G, g = K.lp_rows(A, b, n, K.TLO, T)
    for j in range(N):
        e = [arb(1) if i == j else arb(0) for i in range(N)]
        hi_ = K.certified_max(G, g, l, u, e); lo_ = K.certified_max(G, g, l, u, [-x for x in e])
        if hi_ is not None: u[j] = hi_.min(u[j])
        if lo_ is not None: l[j] = (-lo_).max(l[j])
    U = [K.iv(P[j][0], P[j][1]) for j in range(N + 1)] + [K.iv(l[j], u[j]) for j in range(N)]
    st, KB = K.krawczyk(U, n, iters=12)
    if st == 'empty':
        return dict(status=st)
    t = sum(KB[N + 1:], arb(0))
    PK = K.domain_contract([[KB[j].lower(), KB[j].upper()] for j in range(N + 1)], N)
    c_lo, c_hi = K.C_range(PK, N)
    bound = (-t.lower()).exp().upper() * c_hi
    return dict(status=st, t=[float(t.lower()), float(t.upper())], C=[float(c_lo), float(c_hi)],
                an_upper=float(bound.upper()), margin=float((K.TWO_OVER_E - bound).lower()),
                theta_q_width=max(float(KB[j].rad()) * 2 for j in range(N + 1)))

if __name__ == '__main__':
    n = int(sys.argv[1]); d = json.load(open(sys.argv[2]))
    for f in sorted(d['found'], key=lambda f: -f['an']):
        c = f['cell']
        ok, nodes, reasons, md = K.check_cell(n, c['lo'], c['hi'], None, None, limit=3000)
        vm = value_margin(n, c['lo'], c['hi'])
        print(json.dumps(dict(an_kkt=round(f['an'], 9), depth=f['depth'], prod_reason=f['reason'], indep_excluded=ok,
                              nodes=nodes, reasons=reasons, krawczyk=vm)))
