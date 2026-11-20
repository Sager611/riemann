#!/usr/bin/env python3
"""mlimaudit falsification checks in mpmath (independent of Arb/python-flint):
 (a) 50 extra zeros n in 14401..20311: gamma and |c| = 1/|rho zeta'(rho)| vs the Arb upper bounds in ext_zeros_128.txt;
 (b) zero count: gamma_20311 < 18292.5 < gamma_20312;
 (c) phases theta/pi mod 2 for picked Weisz zeros at mp.prec = 33300 vs lean_dump.json;
 (d) kernel values k(gamma) by a direct breakpoint sweep of int g(x) g(x - t) dx vs lean_dump.json;
 (e) c = 1/(rho zeta'(rho)) for zeros 1, 2, 72, 14400 from Weisz's gamma vs his c strings.
Usage: python3 mp_falsify.py DATADIR"""
import sys, os, json, time
sys.set_int_max_str_digits(0)
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from npzstream import NpzRows
t0 = time.time()
D = sys.argv[1]
W = os.path.join(D, 'wbeni95-mertens-function-bounds-46f5314'); NPZ = os.path.join(D, 'zeta_zeros_14400.npz')
dump = json.load(open(os.path.join(HERE, 'lean_dump.json')))
out = {}
def pm(x):
    return mp.mpf(x.strip().strip('[]').split('+/-')[0].strip())
SKIPAB = '--skipab' in sys.argv
if SKIPAB:
    out = json.load(open(os.path.join(HERE, 'mp_falsify.json')))

# (a) extra zeros
mp.mp.dps = 30
ext = {}
if SKIPAB: ns_skip = True

for line in open(os.path.join(HERE, 'ext_zeros_128.txt')):
    if line.startswith('#'): continue
    import re
    bl = re.findall(r'\[([^ \]]+)(?: \+/- ([^\]]+))?\]', line)
    (gm, _), (cm, cr) = bl[0], bl[1]
    ext[int(line.split()[0])] = (mp.mpf(gm), mp.mpf(cm) + mp.mpf(cr or 0))
ns = [] if SKIPAB else [14401 + 120 * k for k in range(49)] + [20311]
worst_ratio, worst_dg = mp.mpf(0), mp.mpf(0)
for n in ns:
    rho = mp.zetazero(n)
    c_abs = 1 / (abs(rho) * abs(mp.zeta(rho, derivative=1)))
    g_arb, cu = ext[n]
    worst_dg = max(worst_dg, abs(rho.imag - g_arb))
    assert c_abs <= cu * (1 + mp.mpf('1e-20')), ('|c| above Arb bound', n)
    worst_ratio = max(worst_ratio, cu / c_abs - 1)
if not SKIPAB:
    out['a_ext_zeros_checked'] = len(ns)
    out['a_max_|gamma_mp - gamma_arb|'] = mp.nstr(worst_dg, 3)
    out['a_max_(cup/|c|_mp - 1)'] = mp.nstr(worst_ratio, 3)
print('(a) done %.0f s' % (time.time() - t0), flush=True)

# (b) count bracket
if not SKIPAB:
    g1, g2 = mp.zetazero(20311).imag, mp.zetazero(20312).imag
    out['b_gamma_20311'] = mp.nstr(g1, 15); out['b_gamma_20312'] = mp.nstr(g2, 15)
    out['b_bracket_ok'] = bool(g1 < mp.mpf('18292.5') < g2)
print('(b) done %.0f s' % (time.time() - t0), flush=True)
json.dump(out, open(os.path.join(HERE, 'mp_falsify.json'), 'w'), indent=1)

# collect picked Weisz rows (1-based keys of lean_dump.json) without loading the arrays
picks = sorted(int(k) for k in dump)
rows = {}
streams = [NpzRows(NPZ, nm) for nm in ('gamma', 'c_re', 'c_im')]
for i, (gs, cr, ci) in enumerate(zip(*streams)):
    if i + 1 in picks: rows[i + 1] = (gs, cr, ci)
Ys = {w: open(os.path.join(W, 'points', 'y_%s.txt' % w)).read().strip() for w in ('plus', 'minus')}

# (c) phases at 33300 bits
mp.mp.prec = 33300
twopi = 2 * mp.pi
yv = {w: mp.mpf(Ys[w]) for w in Ys}
worst_phase = mp.mpf(0); thetas = {}
for n in picks:
    g = mp.mpf(rows[n][0])
    for w in yv:
        f = g * yv[w] / twopi
        a = 2 * (f - mp.floor(f))
        thetas[(n, w)] = a
        worst_phase = max(worst_phase, abs(a - pm(dump[str(n)][w]['theta_over_pi_mod2'])))
out['c_phase_points'] = len(thetas)
out['c_max_|theta/pi_mp - theta/pi_arb|'] = mp.nstr(worst_phase, 3)
print('(c) done %.0f s' % (time.time() - t0), flush=True)

# (d) kernel values by a breakpoint sweep of the autocorrelation, and the terms 2 k Re(c e^{i theta})
mp.mp.dps = 50
def load_b(w):
    return [int(l) for l in open(os.path.join(HERE, '..', 'MLimsup', 'kernel_%s.txt' % w)) if l.strip() and not l.startswith('#')]
h = mp.mpf(1355) / 100
def k_sweep(b, t):
    M = len(b); T1 = M * h; L = -T1 / 2
    P = sorted([L + j * h for j in range(M + 1)] + [L + j * h + t for j in range(M + 1)])
    acc = mp.mpf(0)
    for p_, q_ in zip(P, P[1:]):
        if q_ - p_ <= 0: continue
        m = (p_ + q_) / 2
        j = int(mp.floor((m - L) / h)); l = int(mp.floor((m - t - L) / h))
        if 0 <= j < M and 0 <= l < M: acc += (q_ - p_) * b[j] * b[l]
    return acc / (h * sum(x * x for x in b))
worst_k, worst_term = mp.mpf(0), mp.mpf(0)
for w in ('plus', 'minus'):
    b = load_b(w)
    for n in picks:
        g = mp.mpf(rows[n][0][:45])
        kv = k_sweep(b, g)
        worst_k = max(worst_k, abs(kv - pm(dump[str(n)][w]['k'])))
        a = thetas[(n, w)]
        term = 2 * kv * (mp.mpf(rows[n][1]) * mp.cospi(a) - mp.mpf(rows[n][2]) * mp.sinpi(a))
        worst_term = max(worst_term, abs(term - pm(dump[str(n)][w]['term_k'])))
out['d_max_|k_sweep - k_arb|'] = mp.nstr(worst_k, 3)
out['d_max_|term_mp - term_arb|'] = mp.nstr(worst_term, 3)
print('(d) done %.0f s' % (time.time() - t0), flush=True)

# (e) c from Weisz's gamma vs his c strings
mp.mp.dps = 75
worst_c = mp.mpf(0)
for n in (1, 2, 72, 14400):
    if n not in rows: continue
    rho = mp.mpc(mp.mpf('0.5'), mp.mpf(rows[n][0][:90]))
    c = 1 / (rho * mp.zeta(rho, derivative=1))
    worst_c = max(worst_c, abs(c.real - mp.mpf(rows[n][1])), abs(c.imag - mp.mpf(rows[n][2])))
out['e_max_|c_mp - c_weisz|'] = mp.nstr(worst_c, 3)
print('(e) done %.0f s' % (time.time() - t0), flush=True)
json.dump(out, open(os.path.join(HERE, 'mp_falsify.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
