
# ---- Good lemma for l = 0, small saddle radius r in [r_lo, r_hi] subset [0, 0.15], parametrized by n >= nmin ----
import fh2lib as L
from fh2lib import arb, acb, ab, ball, polar, up, S, Psie, DkPsie, DkP, fvar, PI, sup_polar, lg_box
from good import lo, hi, ck
def disc0_lg(nu, d0=arb('0.3')):
    # bounds for |h_I|, |h_K| on the closed disc |z| <= d0 (max principle on |z| = d0)
    lgR = lg_box(d0, d0, 0, PI / 2, nu)                 # right half circle (conj-symmetric)
    if lgR is None: return None
    eI, eK, lI, lK = lgR
    # left half: h_I even; 1+h_K(z) = 1+h_K(-z) - kappa e^{2 nu eta(-z)}(1+h_I(-z)), kappa<=1,
    # Re eta(w) <= eta(d0) on |w| = d0, Re w >= 0  (Re eta = Re Psie + log|w| <= Psie(d0) + log d0)
    etamax = Psie(acb(d0)).real + d0.log()
    eK2 = eK + (2 * nu * etamax).exp() * (1 + eI)
    if not (eK2 < 1): return None
    lK2 = arb((-(1 - eK2).log()).upper())
    return eI, eK2, lI, lK2
def g3(z): s = S(z); return (2 + z * z) / s**3
def g4(z): s = S(z); return 1 / s + 3 / s**5
def D4P_f(z):
    u = 1 / (1 + z * z); w = z * z * u
    return 4 * u * w * (1 - 6 * w + 6 * w * w)
def good_bound_l0small(r_lo, r_hi, nmin, nu0, phi1=arb('1.0'), yJ=arb('0.85'), Ktail=40, verbose=False):
    r_lo, r_hi, nu0, nmin = ab(r_lo), ab(r_hi), ab(nu0), ab(nmin)
    fail = lambda msg: (print('FAIL', msg) if verbose else None)
    tt = r_hi + S(acb(r_hi)).real - 1
    n0 = nmin
    if r_lo > 0:
        n_alt = lo(nu0 * (r_lo + S(acb(r_lo)).real - 1))
        if n_alt > n0: n0 = n_alt
    nu_eff = lo(n0 / tt)
    if nu0 > nu_eff: nu_eff = nu0
    d0 = arb('0.3')
    lg = disc0_lg(nu_eff, d0)
    if lg is None: return fail('disc lg')
    eI, eK, lI, lK = lg
    lam = lK                                 # l = 0: log W = -log(1+h_K)
    rmax0 = r_hi * arb('1.01')
    rho = rmax0 / (d0 - rmax0)
    c1, c2, c3, c4 = ck(1, rho), ck(2, rho), ck(3, rho), ck(4, rho)
    kap_lo0 = 1 / (1 + r_hi / 2)
    # exact critical point r_ex: |log(r_ex/r*)| <= drel.  In u = log rho: F(u) = Df, |F(u*)| <= rmax0^2/2 + lam c1,
    # F'(u) = D^2 f >= nu rho - lam c2 >= 0.99 nu r* - lam c2 >= 0.99 n kap_lo0 - lam c2  for rho >= 0.99 r*
    beta_lo_guess = n0 * kap_lo0 * arb('0.99') - lam * c2
    if not (beta_lo_guess > 0): return fail('beta guess')
    drel = (rmax0**2 / 2 + lam * c1) / beta_lo_guess
    if not (drel < arb('0.009')): return fail('drel')      # e^{-0.009} > 0.99 and e^{0.009} < 1.01 (Cauchy disc radii)
    drel_up = drel.exp() - 1                                # r_ex/r* in [1 - drel, 1 + drel_up]
    kap_lo = (1 - drel) * kap_lo0; kap_hi = 1 + drel_up
    rmax = r_hi * (1 + drel_up)
    bl = n0 * kap_lo - lam * c2                                  # beta >= bl (n = n0); scales ~ n
    bh_coef = kap_hi * (1 + rmax)                                # beta <= n*bh_coef + rmax^2 + lam c2
    bh = n0 * bh_coef + rmax**2 + lam * c2
    G3 = sup_polar(g3, r_lo * (1 - drel), rmax, 0, phi1) if r_lo > 0 else sup_polar(g3, arb(0), rmax, 0, phi1)
    G4 = sup_polar(g4, r_lo * (1 - drel), rmax, 0, phi1) if r_lo > 0 else sup_polar(g4, arb(0), rmax, 0, phi1)
    P3 = sup_polar(lambda z: DkP(z, 3), arb(0), rmax, 0, phi1)
    P4 = sup_polar(D4P_f, arb(0), rmax, 0, phi1)
    M3 = n0 * kap_hi * (1 + rmax * G3) + P3 + lam * c3
    M4 = n0 * kap_hi * (1 + rmax * G4) + P4 + lam * c4
    blm = bl - M4 * phi1**2 / 12
    if not (blm > 0): return fail('blm')
    # monotonicity of E in n (n >= n0): blm = A3 n - B3 needs A3 > 0 (B3 >= 0 automatically)
    if not (kap_lo - kap_hi * (1 + rmax * G4) * phi1**2 / 12 > 0): return fail('blm slope')
    T = 2 * (-bl * phi1**2 / 2).exp() * bh.sqrt() / ((2 * PI).sqrt() * bl * phi1)
    E = (M4 / (8 * bl**2) + arb(5) / 24 * M3**2 * bh.sqrt() / blm**arb('3.5') + T) * (bh / bl).sqrt()
    # arc tail: deficit <= n kap_lo (cos phi - 1);  |W ratio| <= e^{2 lam}; prefactor ratio <= 1
    tail = arb(0)
    phis = [phi1 + (PI / 2 - phi1) * k / Ktail for k in range(Ktail + 1)]
    for k in range(Ktail):
        f0, f1 = phis[k], phis[k + 1]
        tail += (f1 - f0) * (n0 * kap_lo * (f0.cos() - 1) + 2 * lam).exp()
    tail = tail / PI
    # imaginary segments [i r_ex, i yJ]: (1/pi) exp(-n kap_lo + lam_seg + lam) / n
    lgs = lg_box(arb(0), yJ, PI / 2, PI / 2, nu_eff)
    if lgs is None: return fail('lgs')
    lseg = lgs[3]
    seg = (-n0 * kap_lo + lseg + lam).exp() / (PI * n0)
    # monotonicity in n of sqrt(beta)*(tail+seg): need n0*kap_lo*(1-cos phi1) > 1
    if not (n0 * kap_lo * (1 - phi1.cos()) > 1): return fail('mono')
    e = E + (tail + seg) * (2 * PI * bh).sqrt()
    # beta/nu <= t(r_hi) kap_hi (1 + rmax) + (rmax^2 + lam c2)/nu_eff   (n = nu t(r*) <= nu t(r_hi))
    bnu = tt * kap_hi * (1 + rmax) + (rmax**2 + lam * c2) / nu_eff
    out = dict(n0=n0, nu_eff=nu_eff, e=hi(e), E=hi(E), tail=hi(tail), seg=hi(seg), drel=drel, lam=lam, bnu=hi(bnu))
    if verbose: print({k: (v.str(4, radius=False) if isinstance(v, arb) else v) for k, v in out.items()})
    return out
