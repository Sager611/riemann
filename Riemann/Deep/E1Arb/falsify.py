# e1arb: falsification of the certificate (cert.py). F1 independent mpmath quadrature; F2 exact cross-term identity and Conrey's
# general form; F3 parameter perturbations; F4 precision; F5 inadmissible inputs refused; F6 float vs Arb; F7 monotonicity; F8 sinh bound.
import json, random, mpmath as mp
from flint import fmpq, fmpq_poly, arb, ctx
from cert import q, buildP, buildQ, checks, kappa, hball, cwang, ie, I01, load
ctx.prec = 512; mp.mp.dps = 90  # monomial Q of degree 49 has coefficients ~1e30: Horner at 40 digits loses ~30
th, nu, R, P, Q, D = load('par_051.json'); c0, k0, X = kappa(nu, R, P, Q)
print("F0 certified kappa(0.51) =", k0.str(20))
# F1: tanh-sinh quadrature in mpmath of C_P, B_P, W2, W2' from the exact coefficients
cf = lambda F: [mp.mpf(int(F[i].p))/int(F[i].q) for i in range(F.degree()+1)]
ev = lambda c, v: mp.polyval(c[::-1], v)
Pc, dPc, Qc, dQc = cf(P), cf(P.derivative()), cf(Q), cf(Q.derivative()); Rm = mp.mpf(int(R.p))/int(R.q)
qs = {'CP': mp.quad(lambda v: ev(dPc, v)**2, [0, 1]), 'BP': mp.quad(lambda v: ev(Pc, v)**2, [0, 1]),
      'W2': mp.quad(lambda v: mp.e**(2*Rm*v)*ev(Qc, v)**2, [0, 0.5, 1]), 'W2p': mp.quad(lambda v: mp.e**(2*Rm*v)*(Rm*ev(Qc, v)+ev(dQc, v))**2, [0, 0.5, 1])}
for kk in qs:
    a = X[kk] if isinstance(X[kk], arb) else arb(X[kk]); print(f"F1 {kk}: quad {mp.nstr(qs[kk], 25)} arb {a.str(25)} rel.diff {mp.nstr(abs(qs[kk]-mp.mpf(a.mid().str(35, radius=False)))/qs[kk], 3)}")
# F2: int e^{2Rv} Q (RQ+Q') = (e^{2R}Q(1)^2 - Q(0)^2)/2 = -1/2 exactly; Conrey's form |P(1)Q(0)|^2 + (1/nu) int int e^{2Rv}|QP' + nu(RQ+Q')P|^2
G1, G0 = ie(Q*(R*Q + Q.derivative()), 2*R); print("F2 cross term: e^{2R}-part", G1, " constant part", G0, " int P P' =", I01(P*P.derivative()))
cg = 1 + (X['CP']/nu)*X['W2'] + nu*X['BP']*X['W2p'] + 2*I01(P*P.derivative())*(arb(G1)*arb(2*R).exp() + arb(G0)); print("F2 Conrey form c =", cg.str(25), " cert c =", c0.str(25))
# F3: perturbations (all admissible); kappa must not exceed the certified optimum by more than rounding
random.seed(1); worst = None
for i in range(24):
    s = lambda e: fmpq(random.randint(-1000, 1000), 1000)*e
    R1 = R*(1 + s(fmpq(1, 1000))) if i % 3 == 0 else R; r1 = q(D['r'])*(1 + s(fmpq(1, 1000))) if i % 3 == 1 else q(D['r'])
    xs = [q(x)*(1 + s(fmpq(1, 10**4))) for x in D['x']] if i % 3 == 2 else [q(x) for x in D['x']]
    P1, Q1 = buildP(r1, 12), buildQ(xs); assert all(checks(th, nu, R1, P1, Q1).values())
    d = kappa(nu, R1, P1, Q1)[1] - k0; worst = d if worst is None or d > worst else worst
print("F3 24 perturbations (R, r +-1e-3 rel; x_j +-1e-4 rel): max kappa - kappa_cert =", worst.str(5))
for nn in ['1324/10000', '13249/100000', '132499/1000000']:
    print("F3 same detector at nu =", nn, ": kappa =", kappa(q(nn), R, P, Q)[1].str(15))
# F4: precision
a1, b1 = ie(Q*Q, 2*R); print("F4 cancellation in W2 = e^{2R}G1 + G0: |G0| ~", arb(b1).str(5), " W2 ~", X['W2'].str(5), " max |monomial coef of Q| ~", f"{float(max(abs(Q[i]) for i in range(Q.degree()+1))):.2e}")
for p in [128, 256, 384, 1024]:
    ctx.prec = p; kp = kappa(nu, R, P, Q)[1]; print(f"F4 prec {p}: kappa = {kp.str(20)}  finite: {kp.is_finite()}  contains-compatible with 512-bit ball: {kp.overlaps(k0) and kp.is_finite()}")
ctx.prec = 512
# F5: inadmissible inputs must be refused
print("F5 nu = 53/400:", {k: v for k, v in checks(th, fmpq(53, 400), R, P, Q).items() if not v})
print("F5 Q + y^2/10^6:", {k: v for k, v in checks(th, nu, R, P, Q + fmpq_poly([0, 0, fmpq(1, 10**6)])).items() if not v})
print("F5 P(1) != 1 (r-Taylor P scaled by 1+1e-9):", {k: v for k, v in checks(th, nu, R, P*(1 + fmpq(1, 10**9)), Q).items() if not v})
# F6, F7: float vs Arb for every theta; dh/dkappa > 0
for t in ['051', '052', '055', '060']:
    tt, n1, R1, P1, Q1, D1 = load(f'par_{t}.json'); k1 = kappa(n1, R1, P1, Q1)[1]; ct = cwang(tt); s = ((1-k1)*(9-k1-8*ct)).sqrt()
    dh = (1 + ((9-k1-8*ct) + (1-k1))/(2*s))/4
    print(f"F6 theta {tt}: float kappa {D1['float_kappa']:.10f} arb {k1.str(12)} diff {float((k1 - arb(D1['float_kappa'])).mid()):.1e}; F7 dh/dkappa = {dh.str(8)}")
# F8: optimiser's closed-form inf over P (sqrt(ab) coth lam) vs the certified polynomial P
a = X['W2']/arb(nu); b = arb(nu)*X['W2p']; lam = (b/a).sqrt(); cs = arb(fmpq(1, 2)) + (a*b).sqrt()/lam.tanh()
print("F8 c(poly P) - c(sinh P) =", (c0 - cs).str(5), " (>= 0: inf over P), lam =", lam.str(12), " r =", D['r'])
# F9: F1 (independent quadrature) and F3-style perturbation for the other certificates
for t in ['052', '055', '060', 'uniform']:
    tt, n1, R1, P1, Q1, D1 = load(f'par_{t}.json'); c1, k1, X1 = kappa(n1, R1, P1, Q1); Rm = mp.mpf(int(R1.p))/int(R1.q)
    Qc, dQc = cf(Q1), cf(Q1.derivative())
    w2 = mp.quad(lambda v: mp.e**(2*Rm*v)*ev(Qc, v)**2, [0, 0.5, 1]); w2p = mp.quad(lambda v: mp.e**(2*Rm*v)*(Rm*ev(Qc, v)+ev(dQc, v))**2, [0, 0.5, 1])
    rd = max(abs(w2 - mp.mpf(X1['W2'].mid().str(40, radius=False)))/w2, abs(w2p - mp.mpf(X1['W2p'].mid().str(40, radius=False)))/w2p)
    random.seed(2); worst = max((kappa(n1, R1*(1 + fmpq(random.choice([-1, 1]), 1000)), P1, Q1)[1] - k1).mid() for i in range(4))
    print(f"F9 {t}: theta {tt} quad-vs-arb max rel.diff (W2, W2') {mp.nstr(rd, 3)}; R +-1e-3 rel: max kappa - kappa_cert {worst.str(5)}; checks {all(checks(tt, n1, R1, P1, Q1, D1.get('uniform', False)).values())}")
