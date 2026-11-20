# b3pair (copy of B3Nu/code/b3nu_arb.py; unchanged but the file name) — b3nu step (3): Arb certificate of the Lemma C threshold with the frozen SL/EXP-010 detector (P, Q, R = 44/5),
# construction executed verbatim from the read-only copy b3c_arb_copy.py (= Riemann/Deep/B3Cert/code/b3c_arb.py) up to kappa().
# Moment range (Lemma C, pair (k,l) = A(k',l')): nu < nu1(theta) = [theta(1+2l) - (k+l)]/(2(1+l-k)); edge theta_c, nu0 = nu1(theta_c).
# usage: b3p_arb.py <params_exp010.json> <prec> [k l theta_c]  (default 1/30 13/15 24701/50000)
import sys
src = open(__file__.replace('b3p_arb.py', 'b3c_arb_copy.py')).read().split('# admissibility')[0]
prec = int(sys.argv[2]); src = src.replace('ctx.prec = 4000', f'ctx.prec = {prec}')
exec(src)
s2 = arb(2).sqrt(); cth = lambda th: 2 - th/2 - (th/s2).cot()/s2; kst = lambda th: 1 - 2/(2 - cth(th))
h = lambda th, kk: (3 + kk - ((1 - kk)*(9 - kk - 8*cth(th))).sqrt())/4
print("detector: P(0) =", P(fmpq(0)), " P(1) =", P(fmpq(1)), " Qp(0) =", Qp(fmpq(0)), " Qp(y)+Qp(1-y)-1 == 0:", (Qp + Qp(fmpq_poly([1, -1])) - 1) == 0, " deg Qp =", Qp.degree(), " prec =", prec, flush=True)
k3, _ = kappa(fmpq(3, 25), fmpq(177, 20)); print("check E1 kappa1 (nu=3/25, R=177/20) =", k3.str(30), flush=True)
qq = lambda t: fmpq(*map(int, t.split('/'))) if '/' in t else fmpq(int(t))
kp, lp, thc = (qq(a) for a in sys.argv[3:6]) if len(sys.argv) > 5 else (fmpq(1, 30), fmpq(13, 15), fmpq(24701, 50000))
R = fmpq(44, 5); nu0 = (thc*(1 + 2*lp) - (kp + lp))/(2*(1 + lp - kp))
print('pair (k,l) =', kp, lp, '; gain factor 1-l-3k =', 1 - lp - 3*kp, flush=True)
print("theta_c =", thc, "=", float(thc), "; nu0 = nu1(theta_c) =", nu0, "=", float(nu0), "; (3 theta_c - 1)/4 =", (3*thc - 1)/4, "; nu0 - (3theta_c-1)/4 =", float(nu0 - (3*thc - 1)/4),
      "== (1-theta_c)(1-l-3k)/(4(1+l-k)):", nu0 - (3*thc - 1)/4 == (1 - thc)*(1 - lp - 3*kp)/(4*(1 + lp - kp)), "; nu0 < theta_c/2:", nu0 < thc/2, flush=True)
k0, c0 = kappa(nu0, R); ks = kst(arb(thc)); hh = h(arb(thc), k0)
print("kappa(nu0, 44/5) =", k0.str(30), "; k*(theta_c) =", ks.str(30), "; kappa - k* =", (k0 - ks).str(20), "(>0:", (k0 - ks) > 0, "); h(theta_c; kappa) =", hh.str(25), "(>0:", hh > 0, ")", flush=True)
kb, _ = kappa((3*thc - 1)/4, R)
print("Cor B' range alone at theta_c: kappa((3theta_c-1)/4) - k*(theta_c) =", (kb - ks).str(15), "(<0:", (kb - ks) < 0, ") -> theta_c is beyond B3Cert's moment edge", flush=True)
print("theta_c < 371/750:", thc < fmpq(371, 750), "; 371/750 - theta_c =", float(fmpq(371, 750) - thc), flush=True)
