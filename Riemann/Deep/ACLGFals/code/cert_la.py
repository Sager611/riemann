# cert_la.py: certified finite-lambda instances of the lower half of Hypothesis LA (j = 1) and of the gap hypothesis, from
# Chekmachev's certified eps_c >= 0.97 Theta_c (Zenodo 23848597, Thm 1.1(b); audited sound in ../../ChekAudit) and weilsimple's
# certified L_odd, L_2 (../../WeilSimple/REPORT.md), against our exact rungs (logs/rungs_cert.json). Theta_c as in ChekAudit/evidence/theta.py.
import json, os, mpmath as mp
mp.mp.dps = 60
here = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(here, '..', 'logs', 'rungs_cert.json')))
def Phi(u):
    return 2*mp.fsum((2*mp.pi**2*n**4*mp.e**(9*u/2) - 3*mp.pi*n**2*mp.e**(5*u/2))*mp.e**(-mp.pi*n**2*mp.e**(2*u)) for n in range(1, 40))
def Theta(c):
    c = mp.mpf(c); a = mp.log(c)/2; kap = 2*mp.pi*c - mp.mpf(9)/2
    return (2*mp.pi*Phi(a)/kap)**4
ws = {}
for l in open(os.path.join(here, '..', '..', 'WeilSimple', 'REPORT.md')):
    p = [x.strip() for x in l.split('|')]
    if len(p) > 6 and p[1].isdigit() and p[2].replace('.', '').isdigit() and 'e' in p[3]: ws[int(p[1])] = tuple(float(x) for x in p[3:6])
rows = []; mins = {}
for c in sorted(int(k) for k in R):
    d4, d6, d8 = (mp.mpf(R[str(c)][n]) for n in ('4', '6', '8'))
    lo = 0.97*Theta(c)/d4
    row = {'c': c, 'cert_eps1_over_d4_lower': float(lo)}
    if c in ws:
        U, Lo, L2 = ws[c]; row.update({'U_over_d4': U/float(d4), 'Lodd_over_d6': Lo/float(d6), 'L2_over_d8': L2/float(d8)})
    rows.append(row)
print(json.dumps(rows))
print('min certified eps1/delta4 lower ratio over c=2..31,45,47:', min(r['cert_eps1_over_d4_lower'] for r in rows),
      'at c =', min(rows, key=lambda r: r['cert_eps1_over_d4_lower'])['c'])
