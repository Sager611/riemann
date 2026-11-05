
# ---- Good lemma (general case, r >= 0.15): rigorous lower bound for J_l(N) ----
# returns relative error e (valid for all nu >= nu0) s.t.
#   J_l(N) >= V * (1+r0^2)^{1/4} * exp(-lamW0) * (1 - e) / sqrt(2 pi beta_hi),   V = exp(nu A_a(r*) - N log r*)
import fh2lib as L
from fh2lib import arb, acb, ab, ball, polar, up, S, Psie, DkPsie, DkP, fvar, PI, sup_polar, lg_box
def lo(x): return arb(x.lower())
def hi(x): return arb(x.upper())
def ck(k, rho):
    if k == 1: return rho
    if k == 2: return rho + 2 * rho**2
    if k == 3: return rho + 6 * rho**2 + 6 * rho**3
    if k == 4: return rho + 14 * rho**2 + 36 * rho**3 + 24 * rho**4
def lamW_per_a(lg, al_hi):
    # lamW = l lamI + (l+1) lamK <= a * ( (lamI+lamK)/2 + alpha |lamK-lamI|/2 )
    eI, eK, lI, lK = lg
    return (lI + lK) / 2 + abs(lK - lI) / 2 * ab(al_hi)
def Jint_lower(rmin, rmax_, f0, f1, K=8):
    # lower bound of int_{f0}^{f1} sin(2psi)/(2 Re sqrt(1+rho^2 e^{2i psi})) dpsi  for rho in [rmin,rmax_]
    tot = arb(0)
    for k in range(K):
        a0 = f0 + (f1 - f0) * k / K; a1 = f0 + (f1 - f0) * (k + 1) / K
        b = ball(a0, a1)
        w = L.S(polar(rmin, rmax_, a0, a1))
        v = (2 * b).sin() / (2 * w.real.upper())
        tot += (a1 - a0) * lo(v)
    return tot
def good_bound_retry(al_lo, al_hi, r_lo, r_hi, nu0, yJ, phi1, Ktail, verbose, Dg_new):
    return good_bound(al_lo, al_hi, r_lo, r_hi, nu0, yJ, phi1, Ktail, verbose, Dg_force=Dg_new)
def good_bound(al_lo, al_hi, r_lo, r_hi, nu0, yJ=arb('0.85'), phi1=None, Ktail=48, verbose=False, nmin=None, Dg_force=None):
    r_lo, r_hi, nu0 = ab(r_lo), ab(r_hi), ab(nu0)
    if nmin is not None and ab(al_hi) == 1:
        # l = 0 : n = nu * t(r*) >= nmin  =>  nu >= nmin / t(r_hi)
        tt = r_hi + L.S(acb(r_hi)).real - 1
        nu_alt = lo(ab(nmin) / tt)
        if nu_alt > nu0: nu0 = nu_alt
    al_lo, al_hi = ab(al_lo), ab(al_hi)
    al = ball(al_lo, al_hi)
    a_lo = 1 / al_hi
    fail = lambda msg: (print('FAIL', msg) if verbose else None)
    Dg = arb('0.02') * r_lo / (1 + r_lo)
    if Dg_force is not None:
        if Dg_force > arb('0.2') * r_lo: return (print('FAIL Dg too large') if verbose else None)
        Dg = Dg_force
    a_lo_ = 1 / al_hi
    if phi1 is None:
        bguess = nu0 * a_lo_ * (al_lo * r_lo + r_lo * r_lo / L.S(acb(r_lo)).real)
        phi1 = arb(3.2) / bguess.sqrt()
        if phi1 < arb('0.3'): phi1 = arb('0.3')
        if phi1 > arb('1.1'): phi1 = arb('1.1')
        phi1 = arb(float(phi1.mid()))
    disc0 = bool(r_hi < arb('0.15'))
    if disc0:
        import goodsmall
        lgd = goodsmall.disc0_lg(nu0, arb('0.3'))
        if lgd is None: return fail('disc0')
        lg1 = lgd
        delta = arb('0.3') - (r_hi + Dg)
    else:
        delta = arb('0.3') * r_lo if r_lo < arb(2) / 3 else arb('0.2')
        th = phi1 + (delta / (r_lo - Dg)).asin()
        lg1 = lg_box(r_lo - delta - Dg, r_hi + delta + Dg, 0, th, nu0)
    if lg1 is None: return fail('lg1')
    lwa = lamW_per_a(lg1, al_hi)                 # lamW <= a*lwa on the Cauchy region
    rho = (r_hi + Dg) / delta
    c1, c2, c3, c4 = ck(1, rho), ck(2, rho), ck(3, rho), ck(4, rho)
    Na0 = nu0 * a_lo
    def Ba_over_a(r): return al * r + r * r / L.S(acb(r)).real
    Bmin = lo(Ba_over_a(r_lo - Dg))
    den = nu0 * Bmin - lwa * c2
    if not (den > 0): return fail('den')
    Delta = (arb('0.5') / a_lo + lwa * c1) * (r_hi + Dg) / den
    if not (Delta < Dg):
        # enlarge the guess and redo (Cauchy region must contain the larger shifts)
        return good_bound_retry(al_lo, al_hi, r_lo, r_hi, nu0, yJ, phi1, Ktail, verbose, Delta * arb('1.4'))
    if not (Delta < Dg): return fail('Delta %s %s' % (Delta.str(3), Dg.str(3)))
    D = Dg
    rmin = r_lo - D; rmax_ = r_hi + D
    bl = lo(Ba_over_a(rmin)) - lwa * c2 / nu0
    bh = hi(Ba_over_a(rmax_)) + (arb('0.25') / a_lo + lwa * c2) / nu0
    m3 = hi(al * sup_polar(lambda z: z, rmin, rmax_, 0, phi1) + sup_polar(lambda z: DkPsie(z, 3), rmin, rmax_, 0, phi1)) \
         + (sup_polar(lambda z: DkP(z, 3), rmin, rmax_, 0, phi1) / a_lo + lwa * c3) / nu0
    m4 = hi(al * sup_polar(lambda z: z, rmin, rmax_, 0, phi1) + sup_polar(lambda z: DkPsie(z, 4), rmin, rmax_, 0, phi1)) \
         + (sup_polar(lambda z: DkP(z, 4), rmin, rmax_, 0, phi1) / a_lo + lwa * c4) / nu0
    blm = bl - m4 * phi1**2 / 12
    if not (blm > 0): return fail('blm')
    Na = Na0
    E = m4 / (8 * Na * bl**2) + arb(5) / 24 * m3**2 * bh.sqrt() / (Na * blm**arb('3.5')) \
        + 2 * (-(Na * bl) * phi1**2 / 2).exp() * bh.sqrt() / ((2 * PI).sqrt() * Na.sqrt() * bl * phi1)
    E = E * (bh / bl).sqrt()
    # contour type: 'up' (r_ex <= yJ: half circle + imaginary segments up to i yJ),
    # 'hor' (r_ex >= yJ: arc to height yJ + horizontal segment), or both (straddling box)
    do_up = bool(rmin < yJ); do_hor = bool(rmax_ > yJ)
    if do_hor and not (rmax_ < arb('1.0') or rmin > yJ):
        pass
    phic = PI / 2 if do_up else (yJ / rmin).asin()
    # the Gaussian arc |phi| <= phi1 must lie on the contour: phi1 <= asin(yJ / r_ex) whenever r_ex > yJ
    if do_hor and not (phi1 < (yJ / rmax_).asin()): return fail('phi1 > phic')
    if not (phi1 < phic): return fail('phi1 >= phic')
    custom = do_hor
    lg0 = lg_box(rmin, rmax_, 0, arb(0), nu0)          # real points: lamW(r)
    lw0 = lamW_per_a(lg0, al_hi)
    # arc tail phi in [phi1, phic]
    Jint = Jint_lower(rmin, rmax_, arb(0), phi1)
    tail = arb(0); worst = None
    phis = [phi1 + (phic - phi1) * k / Ktail for k in range(Ktail + 1)]
    for k in range(Ktail):
        f0, f1 = phis[k], phis[k + 1]
        dpa = hi(al_lo * rmin * (f0.cos() - 1) - rmin**2 * Jint)   # deficit per (nu a), sup over [f0,f1]
        lgk = lg1 if disc0 else lg_box(rmin, rmax_, f0, f1, nu0)
        if lgk is None: return fail('lgk')
        lwk = lamW_per_a(lgk, al_hi)
        # require (lwk+lw0) < nu0 * (-dpa): then exp(a(nu dpa + lwk + lw0)) decreasing in a and nu
        if not (lwk + lw0 < nu0 * (-dpa) / 2): return fail('mono tail %d' % k)
        if not (Na0 * (-dpa) - a_lo * (lwk + lw0) > 1): return fail('mono sqrt tail %d' % k)
        tail += (f1 - f0) * (Na * dpa + a_lo * (lwk + lw0)).exp()
        Jint += Jint_lower(rmin, rmax_, f0, f1, K=2)
        worst = dpa if worst is None or dpa > worst else worst
    tail = tail / PI
    seg = arb(0)
    if do_up:
        dpa = hi(al_lo * rmin * (0 - 1) - rmin**2 * Jint)          # deficit at phi = pi/2 (Jint now = J(pi/2) lower)
        lgs = lg_box(rmin, yJ, PI / 2, PI / 2, nu0)
        if lgs is None: return fail('lgs')
        lws = lamW_per_a(lgs, al_hi)
        if not (lws + lw0 < nu0 * (-dpa) / 2): return fail('mono seg')
        if not (Na0 * (-dpa) - a_lo * (lws + lw0) > 1): return fail('mono sqrt seg')
        tmin = lo(rmin + (L.S(acb(rmin)).real - 1) * a_lo)          # t = rho + a(s-1) >= rmin + a_lo (s(rmin)-1)
        seg += (Na * dpa + a_lo * (lws + lw0)).exp() / (PI * nu0 * tmin)
    if do_hor:
        rlo_h = rmin if rmin > yJ else yJ
        xc = (rmax_**2 - yJ**2).sqrt()
        Kx = 40; segh = arb(0)
        rb = ball(rlo_h, rmax_)
        R = rmax_                              # F(rho, z) is increasing in rho (|z| <= rho): evaluate at rmax
        sR = L.S(acb(R)).real
        for k in range(Kx):
            x0 = xc * k / Kx; x1 = xc * (k + 1) / Kx
            z = acb(ball(x0, x1), yJ)
            lz = (abs(z) / R).log()
            b1 = hi(z.real - R - R * lz)
            b2 = hi(Psie(z).real - Psie(acb(R)).real - (sR - 1) * lz)
            F = (al_lo * b1 if b1 < 0 else al_hi * b1) + b2
            if not (F < 0): return fail('F seg %d' % k)
            pre = hi(abs(1 + z * z)**arb('0.25') / (1 + rb * rb)**arb('0.25') * rb / abs(z))
            ang_lo = lo((yJ / x1).atan())
            ang_hi = hi((yJ / x0).atan()) if x0 > 0 else PI / 2
            if ang_hi > PI / 2: ang_hi = PI / 2
            lgk = lg_box(lo(abs(z)), hi(abs(z)), ang_lo, ang_hi, nu0)
            if lgk is None: return fail('lg seg')
            lwk = lamW_per_a(lgk, al_hi)
            if not (lwk + lw0 < nu0 * (-F) / 2): return fail('mono hseg %d' % k)
            if not (Na0 * (-F) - a_lo * (lwk + lw0) > 1): return fail('mono sqrt hseg %d' % k)
            segh += (x1 - x0) * pre * (Na * F + a_lo * (lwk + lw0)).exp() / lo(rb)
        seg += segh / PI
    beta_hi0 = Na * bh
    e = E + (tail + seg) * (2 * PI * beta_hi0).sqrt()
    out = dict(nu_eff=nu0, phi1=phi1, e=hi(e), E=hi(E), tail=hi(tail), seg=hi(seg), bl=bl, bh=bh, lw0=lw0, lwa=lwa, Delta=Delta, custom=custom, worst=worst)
    if verbose: print({k: (v.str(4, radius=False) if isinstance(v, arb) else v) for k, v in out.items()})
    return out
