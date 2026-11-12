#!/usr/bin/env python3
# falsify.py -- r2aaudit: exact (Arb) tests of Gracia Llorente's Prop 2 identity and Lemma 3 (slot lemma) on actual n,
# and E(x) for primes x <= XMAX against the majorant's x-dependence.  One thread; bounded sizes.
# uv run --no-project --with python-flint --with numpy python3 Riemann/Deep/R2AAudit/falsify.py [quick]
import sys, math, random, time, resource, json, os
import numpy as np
from flint import arb, ctx
ctx.prec = 320
QUICK = 'quick' in sys.argv
DO_ID = 'ex' not in sys.argv; DO_EX = 'id' not in sys.argv     # parts: 'id' (identity/slot tests), 'ex' (E(x) scan); default both
HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time(); LOG = []
def say(*a):
    s = ' '.join(str(v) for v in a); print(s, flush=True); LOG.append(s)
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20
LOGN_MAX = 3.0e4 if QUICK else 2.5e5          # every tested n has log n <= LOGN_MAX
SMAX = int(LOGN_MAX*1.2) + 1000
sv = np.ones(SMAX + 1, dtype=bool); sv[:2] = False
for i in range(2, int(SMAX**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
P = [int(p) for p in np.nonzero(sv)[0]]; del sv
NP = len(P)
EG = arb.const_euler()
# exact prefix sums: TH[k] = theta(p_k) = sum_{i<=k} log p_i, SF[k] = sum_{i<=k} F(p_i), F(p) = -log(1-1/p); index 1-based, [0] = 0
TH = [arb(0)]; SF = [arb(0)]; THf = [0.0]
for p in P:
    lp = arb(p).log(); TH.append(TH[-1] + lp); SF.append(SF[-1] - (1 - arb(1)/p).log()); THf.append(THf[-1] + math.log(p))
say('sieve to %d: %d primes, prefix sums built, %.1f s, %.0f MB' % (SMAX, NP, time.time() - T0, rss()))
def Fq(q): return -(1 - arb(1)/q).log()
def Da(q, a): return -(1 - arb(q)**(-(a + 1))).log()          # D_a(q) = -log(1 - q^{-a-1}) > 0
def analyse(fac):
    """fac: dict prime -> exponent >= 1. Returns dict of exact (ball) quantities and checks."""
    qs = sorted(fac); a = [fac[q] for q in qs]; m = len(qs)
    logn = sum((arb(q).log()*e for q, e in zip(qs, a)), arb(0))
    lf = sum(math.log(q)*e for q, e in zip(qs, a))
    assert lf < THf[-1] - 50, 'n too large for the sieve'
    # k = max index with theta(p_k) <= log n  (float bisection, then exact correction)
    lo_, hi_ = 1, NP
    import bisect
    k = bisect.bisect_right(THf, lf) - 1
    def le(j):   # exact test N_j <= n
        if TH[j] <= logn: return True
        if TH[j] > logn: return False
        nn = 1
        for q, e in zip(qs, a): nn *= q**e
        Nj = 1
        for p in P[:j]: Nj *= p
        return Nj <= nn
    while k + 1 <= NP and le(k + 1): k += 1
    while not le(k): k -= 1
    assert le(k) and not le(k + 1), 'k not certified'
    x = P[k - 1]; T = TH[k]; U = logn - T; L = arb(x).log()
    E = SF[k] - EG - T.log().log()
    sFq = sum((Fq(q) for q in qs), arb(0))
    B = SF[k] - sFq
    D = sum((Da(q, e) for q, e in zip(qs, a)), arb(0))
    H = logn.log().log() - T.log().log()
    lhs = sum((((1 - arb(q)**(-(e + 1)))/(1 - arb(1)/q)).log() for q, e in zip(qs, a)), arb(0)) - EG - logn.log().log()
    resid = lhs - (E - B - D - H)
    slot = sum((Da(P[j], a[j]) for j in range(m)), arb(0))  # sum_j D_{a_j}(p_j)
    r = dict(k=k, x=x, m=m, v2=fac.get(2, 0), lhs=lhs, E=E, B=B, D=D, H=H, resid=resid, slot=slot)
    r['m_le_k'] = m <= k
    r['qj_ge_pj'] = all(qs[j] >= P[j] for j in range(m))
    r['ident_ok'] = resid.contains(arb(0)) and resid.rad() < arb('1e-80')
    r['signs_ok'] = (B > -arb('1e-80')) and (D > 0) and (H > -arb('1e-80'))
    r['slot_margin'] = B + D + H - slot                      # Lemma 3: >= 0 (in fact >= sum c(a_j-1)log p_j + unused slots)
    r['slot_ok'] = r['slot_margin'] > -arb('1e-80')
    r['robin'] = lhs < 0                                    # Robin at n  <=>  lhs < 0
    r['robin_eq'] = (B + D + H > E) if (lhs < 0) else not (B + D + H > E)
    # GL's c = 0.998/(xL) and the chain B+D+cU >= sum_j [D_{a_j}(p_j) + c(a_j-1)log p_j] + sum_{m<i<=k}(F(p_i) - c log p_i)
    c = arb('0.998')/(arb(x)*L)
    r['c_hyp'] = (T + U)*(T + U).log() <= arb(x)*L/arb('0.998')
    if r['c_hyp']:
        rhs = sum((Da(P[j], a[j]) + c*(a[j] - 1)*arb(P[j]).log() for j in range(m)), arb(0)) + (SF[k] - SF[m]) - c*(TH[k] - TH[m])
        r['chain_ok'] = (B + D + c*U - rhs) > -arb('1e-80') and (H - c*U) > -arb('1e-80')
    else:
        r['chain_ok'] = None
    return r
def factor_small(n, spf):
    f = {}
    while n > 1:
        p = int(spf[n]); f[p] = f.get(p, 0) + 1; n //= p
    return f
def sigma_check(fac):
    """independent: sigma(n) as an exact integer, Robin via arb e^gamma n log log n (only for small n)"""
    n = 1; s = 1
    for q, e in fac.items(): n *= q**e; s *= (q**(e + 1) - 1)//(q - 1)
    rhs = EG.exp()*arb(n)*arb(n).log().log()
    return (arb(s) < rhs), (arb(s) > rhs)
STATS = {}
def record(tag, r, robin_direct=None):
    st = STATS.setdefault(tag, dict(n=0, ident_fail=0, sign_fail=0, slot_fail=0, robin_eq_fail=0, mk_fail=0, chain_fail=0, chain_tested=0,
                                    robin_fail=0, direct_mismatch=0, min_slot_rel=None, max_resid=0.0, slot_conclusion_cases=0, slot_conclusion_fail=0))
    st['n'] += 1
    st['ident_fail'] += (not r['ident_ok']); st['sign_fail'] += (not r['signs_ok']); st['slot_fail'] += (not r['slot_ok'])
    st['robin_eq_fail'] += (not r['robin_eq']); st['mk_fail'] += (not (r['m_le_k'] and r['qj_ge_pj']))
    if r['chain_ok'] is not None: st['chain_tested'] += 1; st['chain_fail'] += (not r['chain_ok'])
    st['robin_fail'] += (not r['robin'])
    if not r['robin']:
        st['robin_false_max_x'] = max(st.get('robin_false_max_x', 0), r['x'])      # x = p_k of a Robin failure (5040 = N_4*24: x <= 7)
    st['max_resid'] = max(st['max_resid'], float(abs(r['resid']).upper().mid()) if hasattr(abs(r['resid']).upper(), 'mid') else float(abs(r['resid']).mid()))
    rel = float((r['slot_margin']/(r['B'] + r['D'] + r['H'])).mid())
    if st['min_slot_rel'] is None or rel < st['min_slot_rel'][0]: st['min_slot_rel'] = (rel, r['x'], r['v2'])
    # the conclusion rule: if E(x) < 2^{-v2-1} (v2 >= 1) then Robin must hold
    if r['v2'] >= 1 and r['E'] < arb(2)**(-(r['v2'] + 1)):
        st['slot_conclusion_cases'] += 1; st['slot_conclusion_fail'] += (not r['robin'])
    if robin_direct is not None:
        holds, fails = robin_direct
        if (holds and not r['robin']) or (fails and r['robin']): st['direct_mismatch'] += 1

# ------------------------------------------------------------------ test sets
random.seed(20261112)
def hr_integers(limit, cap):
    """n with non-increasing exponents on the first primes (Hardy-Ramanujan integers), n <= limit (exact integers)"""
    out = []
    def rec(i, n, emax, fac):
        if len(out) >= cap: return
        if i > 0: out.append(dict(fac))
        p = P[i]; pe = n; e = 0
        while e < emax:
            pe *= p; e += 1
            if pe > limit: break
            fac[p] = e; rec(i + 1, pe, e, fac); del fac[p]
    rec(0, 1, 200, {})
    return out
def ca_number(eps):
    """colossally abundant N(eps): exponent of p = floor(log((p^{1+eps}-1)/(p^eps-1))/log p) - 1"""
    fac = {}
    for p in P:
        a = int(math.floor(math.log((p**(1 + eps) - 1)/(p**eps - 1))/math.log(p) + 1e-12)) - 1
        if a <= 0: break
        fac[p] = a
    return fac
def logn_of(fac): return sum(math.log(q)*e for q, e in fac.items())
def run_set(tag, facs, direct=False, budget=None):
    t = time.time(); cnt = 0
    for fac in facs:
        if logn_of(fac) > LOGN_MAX or logn_of(fac) < math.log(6) - 1e-9: continue   # n >= 6: k >= 2, T = log N_k > 1
        r = analyse(fac)
        rd = sigma_check(fac) if direct else None
        record(tag, r, rd); cnt += 1
        if budget and time.time() - t > budget: say('  %s: time budget hit after %d' % (tag, cnt)); break
    st = STATS[tag]
    say('%s: %d n; identity fails %d (max |resid| %.1e); sign fails %d; m<=k,q_j>=p_j fails %d; slot fails %d (min rel margin %.3e at x=%s, v2=%s); '
        'Robin<=>B+D+H>E fails %d; Robin false at %d n; GL chain (c=0.998/(xL)) tested %d fails %d; E<2^-v2-1 cases %d, Robin fails among them %d; direct-sigma mismatches %d; %.0f s, %.0f MB' % (
        tag, st['n'], st['ident_fail'], st['max_resid'], st['sign_fail'], st['mk_fail'], st['slot_fail'], st['min_slot_rel'][0], st['min_slot_rel'][1], st['min_slot_rel'][2],
        st['robin_eq_fail'], st['robin_fail'], st['chain_tested'], st['chain_fail'], st['slot_conclusion_cases'], st['slot_conclusion_fail'], st['direct_mismatch'], time.time() - t, rss()))
if DO_ID:
  # (1) every n in [6, NALL] (exact sigma(n) cross-check)
  NALL = 3000 if QUICK else 40000
  spf = np.zeros(NALL + 1, dtype=np.int64)
  for p in P:
      if p > NALL: break
      sl = spf[p::p]; sl[sl == 0] = p
  run_set('all n in [6,%d]' % NALL, (factor_small(n, spf) for n in range(6, NALL + 1)), direct=True)
  fails = [n for n in range(6, NALL + 1) if not analyse(factor_small(n, spf))['robin']]
  say('  Robin false at n =', fails)
  # (2) Hardy-Ramanujan integers (all SA numbers among them), n <= 10^HRX
  HRX = 12 if QUICK else 20
  hr = hr_integers(10**HRX, 100000)
  hr.sort(key=logn_of)
  best = -1.0; sa = []
  for f in hr:
      v = sum(math.log((q**(e + 1) - 1)/(q**e*(q - 1))) for q, e in f.items())
      if v > best + 1e-15: best = v; sa.append(f)
  say('HR integers <= 1e%d: %d; superabundant among them: %d' % (HRX, len(hr), len(sa)))
  run_set('superabundant <= 1e%d' % HRX, sa, direct=True)
  run_set('HR integers <= 1e%d' % HRX, hr, direct=(HRX <= 12), budget=240)
  # (3) colossally abundant numbers
  cas = {}
  for u in np.arange(0.05, 4.2 if QUICK else 5.0, 0.002):
      f = ca_number(10**(-u))
      if f: cas[tuple(sorted(f.items()))] = f
  cal = sorted(cas.values(), key=logn_of)
  say('CA numbers: %d (largest prime up to %d)' % (len(cal), max(max(f) for f in cal)))
  run_set('colossally abundant', cal)
  # (4) 2^v x odd part of CA numbers, v = 1..80 (near-extremal odd part, prescribed v_2)
  fam = []
  for f in cal[::max(1, len(cal)//60)]:
      for v in (1, 2, 3, 5, 8, 13, 20, 30, 40, 50, 60, 70, 71, 72, 73, 74, 75, 80):
          g = dict(f); g[2] = v; fam.append(g)
  run_set('2^v * oddpart(CA), v<=80', fam, budget=300)
  # (5) random even n with v_2 in [1, 80]
  rnd = []
  for i in range(400 if QUICK else 4000):
      v = random.randint(1, 80); g = {2: v}
      for _ in range(random.randint(0, 40)):
          q = P[min(NP - 1, int(random.expovariate(1/60.0)) + 1)]
          g[q] = g.get(q, 0) + random.choice([1, 1, 1, 2, 2, 3, 4, 6, 10, 25])
      rnd.append(g)
  run_set('random even, v2<=80', rnd)
  # (6) primorial multiples N_k * 2^a * r
  def pm_gen():
      for i in range(300 if QUICK else 3000):
          k = random.randint(2, 2500 if not QUICK else 300); g = {p: 1 for p in P[:k]}
          g[2] += random.randint(0, 80)
          for _ in range(random.randint(0, 5)):
              q = P[random.randint(0, k + 50)]; g[q] = g.get(q, 0) + random.randint(1, 8)
          yield g
  run_set('primorial multiples', pm_gen())
  for tag, st in STATS.items(): say('  %s: largest x = p_k at a Robin failure: %s (n <= 5040 has x <= 11)' % (tag, st.get('robin_false_max_x')))
  del hr, sa, cal, fam, rnd

if DO_EX:
  # ------------------------------------------------------------------ (7) E(x) for primes x <= XMAX vs the majorant's x-dependence
  from math import log, sqrt, exp
  XMAX = 10**5 if QUICK else 10**7
  sv = np.ones(XMAX + 1, dtype=bool); sv[:2] = False
  for i in range(2, int(XMAX**0.5) + 1):
      if sv[i]: sv[i*i::i] = False
  PX = np.nonzero(sv)[0].astype(np.float64); del sv
  g_ = 0.57721566490153286061; SLOW = 2 + g_ - log(4*math.pi); H0 = 3e12; BB = 0.87338196661320752
  TAU = 2.9594e-12*1.0001
  def E1(z):
      # exponential integral E1(z), z > 0 (series for z < 1, continued fraction otherwise; double precision)
      if z < 1:
          s = -g_ - log(z); t = 1.0; k = 1
          while True:
              t *= -z/k; d = -t/k; s += d; k += 1
              if abs(d) < 1e-18: break
          return s
      b = z + 1; c = 1e300; d = 1/b; h = d
      for i in range(1, 400):
          an = -i*i; b += 2; d = 1/(an*d + b); c = b + an/c; de = c*d; h *= de
          if abs(de - 1) < 1e-16: break
      return h*exp(-z)
  def K_int(a, L):  # int_L^oo e^{-a u}(1/u + 1/u^2) du = (1 - a) E1(aL) + e^{-aL}/L
      return (1 - a)*E1(a*L) + exp(-a*L)/L
  # Kahan prefix sums over primes
  sF = 0.0; cF = 0.0; sT = 0.0; cT = 0.0
  maxEM = (0, 0); maxZ = (-1, 0); minZ = (1, 0); minR = (10, 0); maxR = (-10, 0); negE = 0; nchk = 0
  chk_pts = []; last_dec = None
  for idx in range(len(PX)):
      p = PX[idx]
      y = -math.log1p(-1/p) - cF; t = sF + y; cF = (t - sF) - y; sF = t
      y = log(p) - cT; t = sT + y; cT = (t - sT) - y; sT = t
      if p < 1000: continue
      x = p; L = log(x); f = 1/L + 1/L**2; T = sT
      E = sF - g_ - log(log(T))
      if E <= 0: negE += 1
      d = T/x - 1; G = d/L - log(1 + math.log1p(d)/L)
      cpsi = 1 + 1e-8 + 1.01624*x**(-1/6)*L/log(2)
      eps_lo = 2*(L + 2)/(14.13*L*(L + 1)); eps_hi = 2*(L + 2)/(L*(L + 1)*H0)
      M = f*exp(-L/2)*(2*cpsi + SLOW*(1 + eps_lo)) + f*(1 + eps_hi)*TAU*(x**(BB - 1) + x**(-BB))/2 + log(2*math.pi)/(x*L) + max(G, 0.0)
      if E/M > maxEM[0]: maxEM = (E/M, x)
      R = E*sqrt(x)/f
      if R > maxR[0]: maxR = (R, x)
      if R < minR[0]: minR = (R, x)
      if x >= 1e5:
          main = K_int(0.5, L) + K_int(2/3, L)                     # int_x^oo (y^{1/2} + y^{1/3}) w(y) dy
          Z = (E - G + 1/(2*x*L) - log(2*math.pi)/(x*L) - main)/(f*exp(-L/2))
          if Z > maxZ[0]: maxZ = (Z, x)
          if Z < minZ[0]: minZ = (Z, x)
      if nchk % 60000 == 0: chk_pts.append((int(x), E, R))
      nchk += 1
  say('E(x) at %d primes 1e3 <= x <= %.0e: E <= 0 at %d; max E/M_small = %.4f at x = %d (M_small = paper majorant with small-x constants, actual d)' % (nchk, XMAX, negE, maxEM[0], maxEM[1]))
  say('  R(x) = E(x) sqrt(x)/f(L): min %.4f (x=%d), max %.4f (x=%d); paper form bound 2 c_psi + 0.0462(1+eps)' % (minR[0], minR[1], maxR[0], maxR[1]))
  say('  Z(x) = [E - G + 1/(2xL) - log(2pi)/(xL) - int(y^1/2 + y^1/3)w]/(f x^-1/2), x in [1e5, %.0e]: min %.5f (x=%d), max %.5f (x=%d); |zero sum| bound S_low(1+eps_lo) = %.5f' % (
      XMAX, minZ[0], minZ[1], maxZ[0], maxZ[1], SLOW*(1 + 2*(log(1e5) + 2)/(14.13*log(1e5)*(log(1e5) + 1)))))
  say('  samples (x, E, R):', ' '.join('(%d, %.4e, %.4f)' % c for c in chk_pts))
  # cross-check the float E(x) against the exact prefix at a few primes inside the arb sieve
  for x in [P[i] for i in (200, 2000, min(NP - 1, 20000))]:
      k = P.index(x) + 1; Ex = SF[k] - EG - TH[k].log().log()
      j = int(np.searchsorted(PX, x)); 
      say('  exact E(%d) = %s' % (x, Ex.str(12, radius=False)))
STATS2 = {k: {kk: (vv if not isinstance(vv, tuple) else list(vv)) for kk, vv in v.items()} for k, v in STATS.items()}
say('peak RSS %.0f MB, %.0f s' % (rss(), time.time() - T0))
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
json.dump(dict(stats=STATS2, log=LOG), open(os.path.join(HERE, 'logs', 'falsify%s%s.json' % ('_quick' if QUICK else '', '' if (DO_ID and DO_EX) else ('_id' if DO_ID else '_ex'))), 'w'), indent=1)
