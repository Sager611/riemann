# compare2.py (bothaaudit): Botha's Piece IV records vs our RIGOROUS upper bounds (logs/points.jsonl, point.py):
#  q_odd = Uo/mu_odd >= 1 and q_l2 = U2/mu_l2 >= 1 are consistency tests of his certified lower bounds (a value < 1 refutes his record);
#  hybrid step: min(mu_odd, mu_l2)(L_{i+1}) / U_ours(L_i) > 1 (our U in place of his fr_sweep U); his own step ratio for reference.
import json, glob, os, sys
from flint import arb, ctx
ctx.prec = 64
HERE = os.path.dirname(os.path.abspath(__file__))
G = glob.glob(os.path.join(HERE, '..', '..', '..', '.scratch', 'bothaaudit', 'b', '*', 'data', 'piece4', 'glue'))[0]
U = {}; mu = {}
for f in glob.glob(G + '/upper_*.jsonl'):
    for l in open(f):
        r = json.loads(l)
        if r.get('par', 0) == 0 and (r['L'] not in U or arb(r['lam0_upper']) < U[r['L']]): U[r['L']] = arb(r['lam0_upper'])
for f in glob.glob(G + '/g_*.json'):
    nm = os.path.basename(f)[:-5].split('_'); r = json.load(open(f))
    assert r['result'] == 'PASS' and r['Ldec'] == nm[1]
    mu[(nm[1], nm[2])] = arb(r['certified_lower_bound']).lower()
P = {}
for l in open(os.path.join(HERE, 'logs', 'points.jsonl')):
    d = json.loads(l); P[d['L_botha']] = d
Ls = sorted({k[0] for k in mu}, key=float); Ls = [L for L in Ls if float(L) >= 1.70]
f = lambda x: x.str(4, radius=False) if isinstance(x, arb) else str(x)
print('L | c_q | N | U_ours | U_his | U_ours/U_his | Uo/mu_odd | U2/mu_l2 | hybrid step (to next) | his step (to next)')
allok = True; mins = []
for i, L in enumerate(Ls):
    d = P.get(L)
    if not d: continue
    Uu = arb(d['U_up']); Uo = arb(d['Uo_up']); U2 = arb(d['U2_up'])
    num, den = (d['c'].replace('.', ''), 10**len(d['c'].split('.')[1])) if '.' in d['c'] else (d['c'], 1)
    ctx.prec = 256; assert arb(int(num))/den < (2*arb(L)).exp(), ('c_q not below e^{2L}', L); ctx.prec = 64
    qo = Uo/mu[(L, 'odd')]; q2 = U2/mu[(L, 'l2')]
    hs = ms = '-'
    if i + 1 < len(Ls):
        Ln = Ls[i + 1]; m = min(mu[(Ln, 'odd')], mu[(Ln, 'l2')]); hs = m/Uu; ms = m/U[L]
        ok = bool(m > Uu); allok &= ok; mins.append((float(hs.mid()), L))
    print(' | '.join(f(x) for x in [L, d['c'], d['N'], Uu, U[L], Uu/U[L], qo, q2, hs, ms]))
print('hybrid steps all pass:', allok, '; min hybrid ratio', min(mins) if mins else '-')
# the first step from c = 30 (CCMFMM's endpoint, U(30) = 3.27353749418e-150 rigorous) to the first window above it
L0 = Ls[0] if float(Ls[0]) > 1.7006 else Ls[1]
print('step 30 ->', L0, ':', f(min(mu[(L0, 'odd')], mu[(L0, 'l2')])/arb('3.27353749418e-150')))
