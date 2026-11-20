#!/usr/bin/env python3
"""mlimaudit: check ALL 14400 Weisz zeros at 96 bits with Arb (acb.zeta_zeros, own c = 1/(rho zeta'(rho))).
Per zero n: |gamma_W - gamma_n| (index and ~25 digits of the ordinate), delta_n >= |c_W,mid - c_n| rigorously, and whether
Weisz's c ball meets ours. Output: E_w = 2 sum_n |k_w(gamma_n)| delta_n, a rigorous bound for the change of h_K(y)
when Weisz's c midpoints are replaced by the true c (removes the c part of the condition if E_w << margin).
Usage: python3 lowprec_all.py DATADIR [N]"""
import sys, os, time, json
sys.set_int_max_str_digits(0)
from flint import arb, acb, acb_series, ctx
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from npzstream import NpzRows
t0 = time.time()
NPZ = os.path.join(sys.argv[1], 'zeta_zeros_14400.npz')
N = int(sys.argv[2]) if len(sys.argv) > 2 else 14400

def load_kernel(w):
    hdr, b = {}, []
    for line in open(os.path.join(HERE, '..', 'MLimsup', 'kernel_%s.txt' % w)):
        s = line.strip()
        if not s: continue
        if s.startswith('#'):
            if '=' in s:
                k_, v_ = s[1:].split('=', 1); hdr[k_.strip()] = v_.strip()
            continue
        b.append(int(s))
    M = len(b)
    R = [sum(b[j] * b[j + d] for j in range(M - d)) for d in range(M)] + [0, 0]
    return {'M': M, 'hn': int(hdr['h_num']), 'hd': int(hdr['h_den']), 'R': R}
KER = {w: load_kernel(w) for w in ('plus', 'minus')}
def kval(K, g):
    t = g * K['hd'] / K['hn']; d = int(t.mid().floor().unique_fmpz())
    if d >= K['M']:
        assert t >= K['M']; return arb(0)
    assert t >= d and t <= d + 1
    tau = t - d
    return (K['R'][d] * (1 - tau) + K['R'][d + 1] * tau) / K['R'][0]

ctx.prec = 96
E = {w: arb(0) for w in KER}
Wt = {w: arb(0) for w in KER}   # 2 sum |k| |c|: phase sensitivity weight
max_dg, max_delta, bad_index, bad_c = 0.0, 0.0, 0, 0
streams = [NpzRows(NPZ, nm) for nm in ('gamma', 'c_re', 'c_re_rad', 'c_im', 'c_im_rad')]
rows = zip(*streams)
n = 1
while n <= N:
    cnt = min(250, N + 1 - n)
    for j, zz in enumerate(acb.zeta_zeros(n, cnt)):
        gs, cr, crr, ci, cir = next(rows)
        g = zz.imag
        assert zz.real == arb(0.5)
        dg = float(abs(arb(gs[:45]) - g).upper())
        max_dg = max(max_dg, dg)
        if dg > 1e-15: bad_index += 1
        zp = acb_series([zz, 1], prec=2).zeta().coeffs()[1]
        cA = 1 / (zz * zp)
        delta = abs(acb(arb(cr), arb(ci)) - cA).upper()
        max_delta = max(max_delta, float(delta))
        wre = arb(cr) + arb(0, (arb(crr) + arb(1e-58)).upper()); wim = arb(ci) + arb(0, (arb(cir) + arb(1e-58)).upper())
        if not (wre.overlaps(cA.real) and wim.overlaps(cA.imag)): bad_c += 1
        for w in KER:
            kk_ = abs(kval(KER[w], g))
            E[w] += 2 * kk_ * arb(delta)
            Wt[w] += 2 * kk_ * abs(cA)
    n += cnt
    if (n - 1) % 2000 == 0 or n > N:
        print('  %d zeros, %.0f s' % (n - 1, time.time() - t0), flush=True)
res = {'zeros': N, 'max_|gamma_W - gamma_arb|': max_dg, 'index_mismatches(>1e-15)': bad_index,
       'max_delta_c': max_delta, 'c_ball_disagreements': bad_c, 'E_plus': E['plus'].str(5), 'E_minus': E['minus'].str(5),
       'W_plus': Wt['plus'].str(8), 'W_minus': Wt['minus'].str(8)}
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(HERE, 'lowprec_all.json'), 'w'), indent=1)
print('time %.0f s' % (time.time() - t0))
