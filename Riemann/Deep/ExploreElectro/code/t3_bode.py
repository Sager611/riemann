# t3_bode.py (explore22): Kramers-Kronig / Bode test sensitivity. For F with F(s)=F(1-s), the Bode sensitivity integral
# (Balazard-Saias-Yor 1999 for zeta): (1/2pi) int log|zeta(1/2+it)| dt/(1/4+t^2) = sum_{Re rho > 1/2} log|rho/(1-rho)|.
# Each off-line pair (rho, conj rho) adds 2 log|rho/(1-rho)| ~ (2 beta - 1)/|rho|^2: the KK/minimum-phase defect of the line data.
import mpmath as mp, json, sys
sys.path.insert(0, '.')
import ecommon as E
mp.mp.dps = 30
rows = []
zs = [('dh', b, g) for b, g in E.DH_OFF] + [('mv', 0.65, float(E.GBAR))] + [('hyp T=3e12, beta-1/2=1e-3', 0.501, 3e12)]
for name, b, g in zs:
    r = mp.mpc(b, g)
    term = 2 * mp.log(abs(r / (1 - r)))   # rho and its conjugate
    rows.append({'zero': name, 'beta': b, 'gamma': g, 'bsy_term_pair': mp.nstr(term, 6), 'approx_(2b-1)/|rho|^2': mp.nstr((2 * b - 1) / abs(r) ** 2, 6)})
tot = sum(mp.mpf(x['bsy_term_pair']) for x in rows if x['zero'] == 'dh')
print(json.dumps({'rows': rows, 'dh_known_four_total': mp.nstr(tot, 6)}, indent=0))
