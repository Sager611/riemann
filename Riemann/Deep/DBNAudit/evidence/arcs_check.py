# C2: independent exact check of the comparison arcs of Zenodo 23752103 (Section 4), written from the TeX.
#  all arcs: joins, domain, s<=S_i(t0), 0<r<=2(s-o)-1e-12, ell=q0-(2+r q0)dt>0, q1<q0, dyadic q1 (<=2^96)
#  sampled arcs (+ first/last of every row): o>=O_i(q0) exactly, [6/6] Pade outgoing bound exactly
#  float: ideal endpoint vs stored q1 for every arc; re-integration with own rates 2(S_i(t0)-O_i(q0))
#  final: T, q_T, T+q_T/2 == B*, B* < 0.107964910079, startup inequality.
# usage: python3 arcs_check.py PLAN.json TABLE.json FLOW.json.gz [sample_every]
import json, gzip, sys, math, random, time
from fractions import Fraction as F
t0 = time.time()
plan = {r['row']: r for r in json.load(open(sys.argv[1]))['rows']}
table = {r['row']: r for r in json.load(open(sys.argv[2]))}
flow = json.load(gzip.open(sys.argv[3], 'rt'))
every = int(sys.argv[4]) if len(sys.argv) > 4 else 50
recs = flow['records']
print('records', len(recs), 'load_s', round(time.time() - t0, 1), flush=True)
# rows: functional identity plan vs table
for i in range(67):
    p, t = plan[i], table[i]
    assert p['time'] == t['time'] and p['q_interval'] == t['q_interval'], i
    assert p['nodes'] == t['functional']['nodes'], ('nodes', i)
    pj = p.get('joint_jet_overlay') or {}; tj = t['functional'].get('joint_jet_overlay') or {}
    assert pj.get('alpha') == tj.get('alpha') and pj.get('centres') == tj.get('centres'), ('jets', i)
fields = {i: {k: F(v) for k, v in table[i]['fields'].items()} for i in range(67)}
def S(i, t):
    a, b = F(plan[i]['time'][0]), F(plan[i]['time'][1])
    th = (t - a)/(b - a); f = fields[i]
    return f['main_lo'] - f['native_debit'] - (1 - th)*f['first_hi'] - th*f['last_hi']
def O(i, q):
    r = plan[i]
    s = F(0)
    for nd in r['nodes']:
        h, a, w = F(nd['height']), F(nd['shift']), F(nd['weight'])
        u = h*h + a*a - q
        s += 2*w*h*u/(u*u + 4*a*a*q)
    jj = r.get('joint_jet_overlay')
    if jj:
        s *= F(jj['alpha'])
        for c in jj['centres']:
            p = F(c['p']); den = p*p - q
            # direct: (1/m!) d^m/dp^m [1/(p-y)+1/(p+y)] = (-1)^m [(p-y)^-(m+1)+(p+y)^-(m+1)]; symmetric in y -> rational in q
            # compute S_j = u^j+v^j by Newton with A=u+v=2p/den, B=uv=1/den
            A, B = 2*p/den, 1/den
            Sm = [F(2), A]
            for j in range(1, len(c['b']) + 1): Sm.append(A*Sm[j] - B*Sm[j - 1])
            acc = F(0)
            for m, bm in enumerate(c['b']):
                acc += F(bm)*p**(m - 1)*(-1)**m*Sm[m + 1]
            s += 2*acc
    return s
def Ofloat(i, q):
    r = plan[i]; s = 0.0
    for nd in r['nodes']:
        h, a, w = float(F(nd['height'])), float(F(nd['shift'])), float(F(nd['weight']))
        u = h*h + a*a - q; s += 2*w*h*u/(u*u + 4*a*a*q)
    jj = r.get('joint_jet_overlay')
    if jj:
        s *= float(F(jj['alpha']))
        for c in jj['centres']:
            p = float(F(c['p'])); y = math.sqrt(q)
            acc = 0.0
            for m, bm in enumerate(c['b']):
                acc += float(F(bm))*p**(m - 1)*(-1)**m*((p - y)**(-m - 1) + (p + y)**(-m - 1))
            s += 2*acc
    return s
PC = [F(1), F(1, 2), F(5, 44), F(1, 66), F(1, 792), F(1, 15840), F(1, 665280)]
def P(x):
    v = F(0)
    for c in reversed(PC): v = v*x + c
    return v
eta = F(1, 10**12); dt0 = F(1, 2048000)
qinit = F(562500000001, 10**12)
assert F(flow['initial_q']) == qinit and F(flow['initial_sigma']) == F(7, 8)
t, q = F(0), qinit
cnt = dict(arcs=0, exact_own=0, exact_pade=0)
minmarg = dict(rate=None, src=None, own=None, cap=None, ell=None, pade=None)
def upd(k, v):
    if minmarg[k] is None or v < minmarg[k]: minmarg[k] = v
random.seed(7)
prev_row = 0
max_float_excess = 0.0; worst_float = None
rowfirst = set(); rowlast = {}
for k, rec in enumerate(recs):
    if rec['row'] != prev_row: rowlast[prev_row] = k - 1
    if k == 0 or rec['row'] != recs[k - 1]['row']: rowfirst.add(k)
    prev_row = rec['row']
rowlast[prev_row] = len(recs) - 1
special = rowfirst | set(rowlast.values())
qf = float(qinit); tf = 0.0
qf_own = float(qinit)
for k, rec in enumerate(recs):
    i = rec['row']
    t0r, dt, q0, q1 = F(rec['t']), F(rec['dt']), F(rec['qin']), F(rec['qout'])
    s, o, r, ell = F(rec['source_lower']), F(rec['own_upper']), F(rec['rate']), F(rec['linear_positive_lower'])
    assert t0r == t and q0 == q, ('join', k)
    assert dt == dt0
    a, b = F(plan[i]['time'][0]), F(plan[i]['time'][1]); cap = F(plan[i]['q_interval'][1])
    assert F(plan[i]['q_interval'][0]) == 0
    assert a <= t0r and t0r + dt <= b, ('time domain', k)
    assert 0 < q0 <= cap, ('cap', k)
    upd('cap', cap - q0)
    if k > 0: assert i in (recs[k - 1]['row'], recs[k - 1]['row'] + 1), ('row order', k)
    Si = S(i, t0r)
    assert s <= Si, ('source', k); upd('src', Si - s)
    assert r > 0 and r <= 2*(s - o) - eta, ('rate', k); upd('rate', 2*(s - o) - eta - r)
    assert ell == q0 - (2 + r*q0)*dt and ell > 0, ('ell', k); upd('ell', ell)
    assert q1 < q0, ('decrease', k)
    den = q1.denominator; assert den & (den - 1) == 0 and den <= 2**96, ('dyadic', k)
    if k in special or k % every == 0 or random.random() < 0.002:
        Oq = O(i, q0); assert o >= Oq, ('own', k); upd('own', o - Oq); cnt['exact_own'] += 1
        x = r*dt
        bound = (q0 + 2/r)*P(-x)/P(x) - 2/r
        assert bound <= q1, ('pade', k); upd('pade', q1 - bound); cnt['exact_pade'] += 1
    # float ideal endpoint
    rf = float(r); q0f = float(q0)
    ideal = (q0f + 2/rf)*math.exp(-rf*float(dt)) - 2/rf
    exc = float(q1) - ideal
    if exc < -1e-15: raise SystemExit(('float ideal above stored q1', k, exc))
    if exc > max_float_excess: max_float_excess = exc; worst_float = k
    # own-rate re-integration (float): r_own = 2(S_i(t0) - O_i(q_own))
    r_own = 2*(float(Si) - Ofloat(i, qf_own))
    qf_own = (qf_own + 2/r_own)*math.exp(-r_own*float(dt)) - 2/r_own
    t, q = t0r + dt, q1
    cnt['arcs'] += 1
T = t; qT = q
assert T == F(22101, 204800) and qT == F(7902377963830536747158643, 2**96) == F(flow['terminal_q'])
B = T + qT/2
Bstar = F(427693072077942473661578098491, 3961408125713216879677197516800)
assert B == Bstar == F(flow['endpoint']) and Bstar < F(107964910079, 10**12)
r0 = F(recs[0]['rate'])
start = F(1, 10**12) - r0*qinit*F(1, 10**14)
assert start > F(1, 2*10**12)
print(json.dumps(dict(counts=cnt, min_margins={k: float(v) for k, v in minmarg.items()},
    max_float_rounding_excess=max_float_excess, worst_float_arc=worst_float,
    terminal=dict(T=str(T), qT=float(qT), endpoint=str(B), endpoint_float=float(B)),
    own_rate_float_reintegration_qT=qf_own, own_rate_endpoint=float(T) + qf_own/2,
    startup_margin=float(start), seconds=round(time.time() - t0, 1))), flush=True)
print('ARCS_CHECK_PASS')
