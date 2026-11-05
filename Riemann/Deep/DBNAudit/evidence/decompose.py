# D: gain decomposition for Zenodo 23752103 vs 23272725 (float model of the published comparison scheme).
# Model: arcs of length delta, rate r = 2(S_i(t0) - O_i(q0)) - 1e-12 at the arc start, exact ODE Q'=-2-rQ, stop at q<=1e-4,
# bound = t + q/2.  Counterfactuals remove one ingredient at a time (sources and own pairs from the author's own fields).
# usage: python3 decompose.py PLAN TABLE DIRECT OLD_HISTORY.json.gz
import signal; signal.alarm(1500)  # hard time limit
import json, gzip, sys, math
from fractions import Fraction as F
plan = {r['row']: r for r in json.load(open(sys.argv[1]))['rows']}
table = {r['row']: r for r in json.load(open(sys.argv[2]))}
direct = {r['row']: r for r in json.load(open(sys.argv[3]))['rows']}
old = json.loads(gzip.decompress(open(sys.argv[4], 'rb').read()))
fl = lambda v: float(F(v))
rows = []
for i in range(67):
    p, t, d = plan[i], table[i], direct[i]
    f = {k: fl(v) for k, v in t['fields'].items()}
    hc = d['head_construction']
    if 'raw_head_bounds' in hc:
        Hm, Hp = [fl(x) for x in hc['raw_head_bounds']]
        Cred = fl(hc['composite_credit']) + fl(hc['old_angular_credit']) + fl(hc['additional_clock_credit'])
    else:
        Hm, Hp = [fl(x) for x in hc['raw_cap_sums']]
        Cred = fl(hc['fixed_credit']) + fl(hc['coordinate_credit']) + fl(hc['orbit_credit'])
    jj = p.get('joint_jet_overlay')
    alpha = fl(jj['alpha']) if jj else 1.0
    m_sc = fl(d['scalar_main_lower']); D_sc = fl(d['scalar_native_before_refund']) - fl(d['refund'])
    jet_main = f['main_lo'] - alpha*m_sc; jet_nat = f['native_debit'] - alpha*D_sc
    jet_hm = f['first_hi'] - alpha*(Hm - Cred); jet_hp = f['last_hi'] - alpha*(Hp - Cred)
    nodes = [(fl(n['height']), fl(n['shift']), fl(n['weight'])) for n in p['nodes']]
    cents = [(fl(c['p']), [fl(b) for b in c['b']]) for c in jj['centres']] if jj else []
    rows.append(dict(a=fl(p['time'][0]), b=fl(p['time'][1]), f=f, alpha=alpha, Cred=Cred, gainM8=fl(t['gain_from_M8_source']) if t.get('gain_from_M8_source') not in (None, '') else 0.0,
                     jet=(jet_main, jet_nat, jet_hm, jet_hp), nodes=nodes, cents=cents, cap=fl(p['q_interval'][1]), Hm=Hm, Hp=Hp, msc=m_sc, Dsc=D_sc))
def O_old(R, q):
    s = 0.0
    for h, a, w in R['nodes']:
        u = h*h + a*a - q; s += 2*w*h*u/(u*u + 4*a*a*q)
    return s
def O_jet(R, q):
    y = math.sqrt(q); s = 0.0
    for p, b in R['cents']:
        for m, bm in enumerate(b):
            s += 2*bm*p**(m - 1)*(-1)**m*((p - y)**(-m - 1) + (p + y)**(-m - 1))
    return s
def S_src(R, t, mode):
    f = R['f']; th = (t - R['a'])/(R['b'] - R['a'])
    S = f['main_lo'] - f['native_debit'] - (1 - th)*f['first_hi'] - th*f['last_hi']
    if 'nocredit' in mode: S -= R['alpha']*R['Cred']
    if 'nogainM8' in mode: S -= R['gainM8']
    if 'nojet' in mode and R['cents']:
        jm, jn, jh, jp = R['jet']
        S -= jm - jn - (1 - th)*jh - th*jp
        # rescale scalar part from alpha to 1
        Ssc = (S - 0)  # scalar part currently alpha*(...)
        S = S/R['alpha']
    return S
def O_tot(R, q, mode):
    if 'nojet' in mode and R['cents']: return O_old(R, q)
    return R['alpha']*O_old(R, q) + (O_jet(R, q) if R['cents'] else 0.0)
def run(mode=(), delta=1/2048000, q0=0.562500000001, qstop=1e-4):
    q, i, n = q0, 0, 0
    ends = [round(R['b']/delta) for R in rows]
    while True:
        while i < 66 and n >= ends[i]: i += 1
        R = rows[i]; t = n*delta
        if n + 1 > ends[i]:
            return dict(t=t, q=q, bound=t + q/2, arcs=n, ended='rows')
        r = 2*(S_src(R, t, mode) - O_tot(R, q, mode)) - 1e-12
        q1 = (q + 2/r)*math.exp(-r*delta) - 2/r
        if q1 <= qstop:
            # exact travel time to qstop inside this arc (continuous stop), then bound
            tau = math.log((q + 2/r)/(qstop + 2/r))/r
            return dict(t=t + tau, q=qstop, bound=t + tau + qstop/2, arcs=n + 1, ended='q', rate_last=r)
        q = q1; n += 1
out = {}
out['base_model'] = run()
out['base_model_grid_stop'] = None
out['no_arith_credits'] = run(('nocredit',))
out['no_post_M8_gains'] = run(('nogainM8',))
out['no_jets'] = run(('nojet',))
out['coarse_grid_1_128000'] = run(delta=1/128000)
out['seed_9_16_plus_1e-6'] = run(q0=0.562501)
out['stop_at_1e-8'] = run(qstop=1e-8)
out['all_but_design'] = run(('nocredit', 'nojet'), delta=1/128000, q0=0.562501, qstop=1e-8)
# old chain (23272725 main primary): rate vs q, and time-to-reach profile
recs = old['records']
tq_old = [(fl(r['t_lo']), fl(r['q_in']), fl(r['rate'])) for r in recs]
# new chain profile from model
def profile(mode=()):
    delta = 1/2048000; q, i, n = 0.562500000001, 0, 0; prof = [(0.0, q, None)]
    ends = [round(R['b']/delta) for R in rows]
    while True:
        while i < 66 and n >= ends[i]: i += 1
        R = rows[i]; t = n*delta
        r = 2*(S_src(R, t, mode) - O_tot(R, q, mode)) - 1e-12
        q1 = (q + 2/r)*math.exp(-r*delta) - 2/r
        prof.append((t + delta, q1, r)); q = q1; n += 1
        if q1 <= 1e-4 or n >= ends[66]: return prof
pnew = profile()
def time_at(prof, qlev):
    for k in range(1, len(prof)):
        if prof[k][1] <= qlev:
            (t0, q0, _), (t1, q1, _) = prof[k - 1], prof[k]
            return t0 + (t1 - t0)*(q0 - qlev)/(q0 - q1)
    return None
pold = [(t, q, r) for (t, q, r) in tq_old]
levels = [0.5625, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001, 0.0001]
band = []
for a, b in zip(levels[:-1], levels[1:]):
    to = time_at(pold, b) - (time_at(pold, a) if a < 0.5625 else 0.0)
    tn = time_at(pnew, b) - (time_at(pnew, a) if a < 0.5625 else 0.0)
    # representative rates at band midpoint
    qm = (a + b)/2
    ro = next(r for (t, q, r) in pold if q <= qm); rn = next(r for (t, q, r) in pnew[1:] if q <= qm)
    band.append(dict(q_from=a, q_to=b, dt_old=to, dt_new=tn, saving=to - tn, rate_old_mid=ro, rate_new_mid=rn))
old_end = {'t_at_q=1e-4': time_at(pold, 1e-4), 'steps': len(recs), 'headline': 0.116397726408291}
out['old_vs_new_by_q_band'] = band
out['old_chain'] = old_end
print(json.dumps(out, indent=1))
