
# large_tau.py : the regime tau = n/nu >= TAU0 (nu >= nu0), explicit linear-in-tau bounds
import sys, json, math
import assembly as A
from fh2lib import arb, acb, ab, ball, polar, S, Psie, PI, lg_box
lo, hi = A.lo, A.hi
nu0 = arb(sys.argv[1]); A.load(nu0, sys.argv[2], sys.argv[3])
TAU0 = arb(sys.argv[4]) if len(sys.argv) > 4 else arb(60)
YJ = A.YJ
def eta_r(r): return Psie(acb(r)).real + arb(r).log()
ok = True
def check(cond, msg):
    global ok
    print(('ok   ' if cond else 'FAIL ') + msg, flush=True)
    ok = ok and bool(cond)
# ---- (a) choice of l-hat: r(25) >= 0.5 at tau = TAU0, and |d rho/d l| <= 2 s^2/(s + a rho) <= 0.1055 on [0.5, 0.635], a >= 51
r25 = A.saddle(51, TAU0 - 50)
check(r25 > arb('0.5'), 'r(l=25, tau=%s) = %s >= 0.5' % (TAU0, r25.str(5)))
sm = S(acb(arb('0.635'))).real
step = 2 * sm**2 / (S(acb(arb('0.5'))).real + 51 * arb('0.5'))
check(step < arb('0.135'), 'saddle step bound 2 s^2/(s + a rho) <= %s < 0.635-0.5' % step.str(4))
# ---- (a2) classification facts used below
s16_ = S(acb(arb('1.6'))).real
check(Psie(acb(arb('1.6'))).real + arb('1.6').log() > 0, 'eta(1.6) > 0 (far-arc rate increasing in a)')
check(all(arb('1.6') + (2 * l + 1) * (s16_ - 1) < TAU0 - 2 * l for l in range(5)), 'l <= 4: r_l > 1.6 for tau >= TAU0 (all far)')
rows5 = sorted(A.TAB['>=5'], key=lambda x: x[0])
cov = [x for x in rows5 if x[1] >= A.RG5 and x[0] <= A.RMAX]
contig = bool(cov) and cov[0][0] <= A.RG5 and cov[-1][1] >= A.RMAX and all(v[0] <= u[1] for u, v in zip(cov, cov[1:]))
check(contig and all(x[2] is not None and x[2] < 0.9 for x in cov), 'family l>=5: Good rows cover r in [%.3f, %.1f] with e < 0.9 (those J_l > 0)' % (A.RG5, A.RMAX))
# ---- (b) main term: Phi_hat >= I0 + g_min * tau on r in [0.5, 0.635]
def gfun(r):
    s = S(acb(r)).real
    return (Psie(acb(r)).real - (s - 1) * r.log()) / s
K = 200; gmin = None; I0 = None
for k in range(K):
    rb = ball(arb('0.5') + arb('0.135') * k / K, arb('0.5') + arb('0.135') * (k + 1) / K)
    gv = lo(gfun(rb))
    iv = lo(rb * (1 - rb.log()) + (1 - rb) * gfun(rb))
    gmin = gv if gmin is None or gv < gmin else gmin
    I0 = iv if I0 is None or iv < I0 else I0
print('Phi_hat >= %s + %s * tau' % (I0.str(5), gmin.str(5)))
# good factor: e, bh, lw0 from family table rows covering [0.5, 0.635]; l-hat = (a-1)/2 with a <= tau+1
rows = [x for x in A.TAB['>=5'] if x[1] > 0.5 and x[0] < 0.635]
emax = max(x[2] for x in rows); bhmax = max(x[3] for x in rows)
lam = A.lam_region(0.4, 0.8, 0, 0.001); lw = (lam[0] + lam[1]) / 2 + abs(lam[1] - lam[0]) / 2  # per unit a
print('good rows: e <= %.4f, beta/(nu a) <= %.4f, lamW0 per a <= %s' % (emax, bhmax, lw.str(4)))
# log gf >= log(1-e) - a lw - l/(12 nu0) - (1/2) log(2 pi nu a bh), a <= tau+1, l <= tau/2
def gf_lower(tau):
    a = tau + 1
    return (1 - arb(emax)).log() - a * lw - (tau / 2) / (12 * nu0) - (2 * PI * nu0 * a * arb(bhmax)).log() / 2
# ---- classes: rate(tau) <= C0 + C1 tau  (per nu), prefactor log <= c0 + c1 tau (+ log count)
classes = []
# Bad-far arc (r_l > 1.6 <=> a < (tau-0.6)/s16)
s16 = S(acb(arb('1.6'))).real; e16 = eta_r(arb('1.6'))
C1 = -arb('1.6').log() + e16 / s16
C0 = arb('1.6') - arb('1.6').log() - e16 * arb('0.6') / s16
lamf = A.lam_region(1.598, 1.602, 0, float((YJ / arb('1.6')).asin()) + 1e-6)
classes.append(('far-arc', C0, C1, hi((1 + arb('1.6')**2).log() / 4) , (lamf[0] + lamf[1]), 'count_far'))
# Bad-far horizontal pieces at height yJ, x in [0, xc]
xc = (arb('1.6')**2 - YJ**2).sqrt(); Kx = 60
hc0 = None; hc1 = None; hcc = None; hlam = arb(0)
for k in range(Kx):
    x0 = xc * k / Kx; x1 = xc * (k + 1) / Kx
    z = acb(ball(x0, x1), YJ)
    lz = abs(z).log(); reta = Psie(z).real + lz
    base0 = z.real - lz                        # Re z - (tau+1) log|z| = (Re z - log|z|) - tau log|z|
    # a in [1, (tau-0.6)/s16]: endpoint a=1: base + Re eta ; endpoint a_max: base + (tau-0.6)/s16 Re eta
    for (c0, c1) in [(base0 + reta, -lz), (base0 - arb('0.6') / s16 * reta, -lz + reta / s16)]:
        c0h, c1h = hi(c0), hi(c1)
        # combine as max over tau >= TAU0 of c0 + c1 tau  <= (c0 + c1 TAU0) + max(c1,0)(tau - TAU0)  -> store as line
        hc0 = c0h if hc0 is None or c0h > hc0 else hc0
        hc1 = c1h if hc1 is None or c1h > hc1 else hc1
    cc = hi(abs(1 + z * z).log() / 4 - abs(z).log())
    hcc = cc if hcc is None or cc > hcc else hcc
    angl = float(lo((YJ / x1).atan())); angh = float(hi((YJ / x0).atan())) if x0 > 0 else 1.5708
    lm = A.lam_region(float(lo(abs(z))), float(hi(abs(z))), angl, min(angh, 1.5708))
    hlam = max(hlam, lm[0] + lm[1])
classes.append(('far-hor', hc0, hc1, hcc + (xc / PI).log(), hlam, 'count_far'))
# Bad-start, contour radius rs = 0.116 (r_l < rs <=> t < rs + a(s_s - 1)), a <= tau + 1
rs = arb('0.116'); ss = S(acb(rs)).real; L = -rs.log()
co = Psie(acb(rs)).real + L * (ss - 1)
check(co > 0 and Psie(acb(0, rs)).real + L * (ss - 1) > 0, 'start-class a-coefficients positive (a <= tau + 1 bound valid)')
classes.append(('start-arc', rs * (1 + L) + co, co, hi((1 + rs * rs).log() / 4), None, 'count_start'))
co2 = Psie(acb(0, rs)).real + L * (ss - 1)
classes.append(('start-seg', L * rs + co2, co2, hi((YJ / rs).log().log() - PI.log()), None, 'count_start'))
lamsA = A.lam_region(0.1155, 0.1165, 0, 1.5708); lamsS = A.lam_region(0.116, 0.85, 1.5707, 1.5708)
lam_start = max(lamsA[0] + lamsA[1], lamsS[0] + lamsS[1])
# Neg: rate <= (tau+3) Psie(0)
P0 = Psie(acb(0)).real
lamN = A.lam_region(0.0, 0.85, 1.5707, 1.5708)
reta = Psie(acb(0, YJ)).real + YJ.log(); lnX = 2 * nu0 * reta + lamN[0] + lamN[1]
check(lnX < 0, 'X(i yJ) < 1')
classes.append(('neg', 3 * P0, P0, hi(-(1 - lnX.exp()).log() + YJ.log() - PI.log()), lamN[0] + lamN[1], 'one'))
# I_L: rate <= M_IL + tau*(-log yJ)
MIL = None; cIL = None
for k in range(48):
    f0 = PI / 2 + PI / 2 * k / 48; f1 = PI / 2 + PI / 2 * (k + 1) / 48
    z = polar(YJ, YJ, f0, f1)
    v = hi((z + Psie(z)).real); MIL = v if MIL is None or v > MIL else MIL
    lm = A.lam_region(0.8499, 0.8501, float(PI - f1), float(PI - f0))
    c = hi(abs(1 + z * z).log() / 4 - (1 - (1 - (-lm[1]).exp())).log())
    cIL = c if cIL is None or c > cIL else cIL
classes.append(('IL', MIL, -YJ.log(), cIL, arb(0), 'one'))
# ---- final inequality: for tau >= TAU0,  sum_c count_c * exp(nu (R_c - Phi) + c_c) < 1/2 exp(gf)
def count(kind, tau):
    if kind == 'count_far': return tau / (2 * s16) + 1
    if kind == 'count_start': return tau / 2 + 1
    return arb(1)
tot = arb(0); worstD = None; slopes_ok = True
for (name, C0, C1, cpre, lamc, kind) in classes:
    D0 = I0 + gmin * TAU0 - (C0 + C1 * TAU0)     # Phi_hat - R at TAU0
    D1 = gmin - C1                               # slope in tau
    if lamc is None: lamc = lam_start
    # prefactor exponent: cpre + lamc*(tau+1) (lambda_W <= a * lamc, a <= tau+1) + log count ; minus gf
    expo0 = -nu0 * D0 + cpre + lamc * (TAU0 + 1) + count(kind, TAU0).log() - gf_lower(TAU0)
    # d/dtau of exponent <= -nu0 D1 + lamc + 1/TAU0 + (d/dtau of -gf <= lw + 1/(24 nu0) + 1/(2 TAU0))
    dexp = -nu0 * D1 + lamc + 1 / TAU0 + lw + 1 / (24 * nu0) + 1 / (2 * TAU0)
    good = bool(D1 > 0) and bool(dexp < 0) and bool(nu0 * D0 > 1)
    slopes_ok = slopes_ok and good
    print('  %-10s rate <= %s + %s tau ; gap at tau0 %s, slope %s ; log(term/main) at tau0 <= %s ; d/dtau <= %s %s' % (
        name, hi(C0).str(4), hi(C1).str(4), lo(D0).str(4), lo(D1).str(4), hi(expo0).str(5), hi(dexp).str(4), 'ok' if good else 'FAIL'))
    tot += hi(expo0).exp()
check(slopes_ok, 'all classes: positive gap, gap slope > 0, exponent decreasing in tau (so worst case is tau = TAU0)')
check(tot < arb('0.5'), 'sum of all bad contributions / main term <= %s < 1/2 at tau = TAU0 (all tau >= TAU0 by monotonicity)' % tot.str(3))
print('LARGE-TAU ALL OK' if ok else 'LARGE-TAU FAILURES')
