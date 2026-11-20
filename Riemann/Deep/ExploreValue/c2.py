# C2 (bounded scale): window [T1,T2]; zeta or DH. (1) N0bar, N_tot, N+; (2) #zeros of V_c right of 1/2;
# (3) Levinson V_b = F + (b/L)F' right of 1/2; (4) Cauchy(1/2) law of eta = Im(F'/F)/|chi'/chi| on the line.
import numpy as np, math, sys, time
from fz import Fn, line_eval
kind = sys.argv[1]; T1 = float(sys.argv[2]); T2 = float(sys.argv[3]); dt = float(sys.argv[4])
f = Fn(kind, T2 + 20); n = int(round((T2 - T1)/dt)) + 1; t = T1 + dt*np.arange(n); t0c = time.time()
F, Fp = line_eval(f, 0.5, T1, n, dt, Bk=50); s = 0.5 + 1j*t; ch = f.chilog(s)
Z = np.real(np.exp(1j*f.theta(t))*F); n0 = int(np.sum(np.sign(Z[1:]) != np.sign(Z[:-1])))
print(f"{kind} [{T1},{T2}] dt={dt} M={len(f.coef)} line eval {time.time()-t0c:.0f}s; Im(chi'/chi) max {np.max(np.abs(ch.imag)):.1e}", flush=True)
sg = np.arange(0.5, 3.0, 0.0005)
def hor(T):
    Fh, Fph = f(sg + 1j*T, blk=40); return Fh, Fph, f.chilog(sg + 1j*T)
Fb, Fpb, chb = hor(T1); Ft, Fpt, cht = hor(T2)
nr = int(round((T2 - T1)/0.05)) + 1; tr = T1 + 0.05*np.arange(nr); Fr, Fpr = line_eval(f, 3.0, T1, nr, 0.05); chr_ = f.chilog(3.0 + 1j*tr)
def S_of(Fh):   # (1/pi) arg F(1/2+iT), continuous from sigma=3
    a = np.unwrap(np.angle(Fh[::-1])); return a[-1]/math.pi
Ntot = (f.theta(np.array([T2]))[0] - f.theta(np.array([T1]))[0])/math.pi + S_of(Ft) - S_of(Fb)
npl = (Ntot - n0)/2
print(f"  line zeros N0bar={n0}  all zeros N_tot={Ntot:.3f}  N+={npl:.2f}", flush=True)
def wind(G_b, G_r, G_t, G_l):   # counterclockwise: bottom(sig up), right(t up), top(sig down), left(t down)
    seq = np.concatenate([G_b, G_r, G_t[::-1], G_l[::-1], G_b[:1]]); ph = np.unwrap(np.angle(seq))
    return (ph[-1] - ph[0])/(2*math.pi), np.max(np.abs(np.diff(ph)))
res = []
for c in [2, 0.5, 0.2, -0.2, -0.4]:
    Vc = lambda F_, Fp_, ch_: Fp_ - (0.5 + c)*ch_*F_
    w, st = wind(Vc(Fb, Fpb, chb), Vc(Fr, Fpr, chr_), Vc(Ft, Fpt, cht), Vc(F, Fp, ch)); res.append(f"c={c:+}: {w:.2f} (step {st:.2f})")
print("  #zeros V_c right of 1/2: " + "; ".join(res), flush=True)
Lm = math.log((T1+T2)/2/(2*math.pi)); res = []
for b in [1.03217, 1.5, 3.0, 6.0]:
    Vb = lambda F_, Fp_: F_ + (b/Lm)*Fp_
    w, st = wind(Vb(Fb, Fpb), Vb(Fr, Fpr), Vb(Ft, Fpt), Vb(F, Fp)); res.append(f"b={b}: {w:.2f} (step {st:.2f})")
print(f"  Levinson V_b = F+(b/L)F', L={Lm:.4f}, #zeros right of 1/2: " + "; ".join(res), flush=True)
eta = np.imag(Fp/F)/np.abs(ch.real); rp = np.real(Fp/F)/np.abs(ch.real)
print(f"  check Re(F'/F)/|chi'/chi| on line: mean {rp.mean():.6f} max|+0.5| {np.max(np.abs(rp+0.5)):.1e}")
out = []
for c in [0.05, 0.1, 0.25, 0.5, 1, 2, 4]:
    D = np.mean(np.log(np.abs(c + 1j*eta))) - math.log(c + 0.5); P = np.mean(c/(c*c + eta*eta)) - 1/(c + 0.5)
    out.append(f"c={c}: Dlog={D:+.5f} Dres={P:+.5f}")
print("  Cauchy(1/2) test, Dlog = E log|c+i eta| - log(c+1/2), Dres = E c/(c^2+eta^2) - 1/(c+1/2):\n   " + "\n   ".join(out))
xs = np.sort(eta); cdf = 0.5 + np.arctan(2*xs)/math.pi; emp = (np.arange(n) + 0.5)/n
print(f"  KS distance to Cauchy(0,1/2): {np.max(np.abs(emp - cdf)):.5f}; median {np.median(eta):+.4f}; P(|eta|<1/2) = {np.mean(np.abs(eta) < 0.5):.4f} (Cauchy 0.5)")
# Littlewood on [1/2,3]x[T1,T2]: 2pi*W(h) = int log|h(1/2+it)| - int log|h(3+it)| + int_top arg h - int_bottom arg h,
# arg by continuous variation leftward from sigma=3 (principal there; |h(3+it)/h(3+iT1)-1| small).  W(h) = sum over zeros right of 1/2 of (beta-1/2).
def W_of(hl, hr, hb, ht):
    ab = np.unwrap(np.angle(hb[::-1]))[::-1]; at = np.unwrap(np.angle(ht[::-1]))[::-1]
    at = at - at[-1] + np.angle(ht[-1]); ab = ab - ab[-1] + np.angle(hb[-1])
    return (np.trapezoid(np.log(np.abs(hl)), t) - np.trapezoid(np.log(np.abs(hr)), tr) + np.trapezoid(at, sg) - np.trapezoid(ab, sg))/(2*math.pi)
out = [f"W(F)={W_of(F, Fr, Fb, Ft):.4f}"]
for c in [2, 0.5, 0.2]:
    Vc = lambda F_, Fp_, ch_: Fp_ - (0.5 + c)*ch_*F_
    out.append(f"W(V_c={c})={W_of(Vc(F,Fp,ch), Vc(Fr,Fpr,chr_), Vc(Fb,Fpb,chb), Vc(Ft,Fpt,cht)):.4f}")
print("  Littlewood weighted counts (sum of beta-1/2 right of the line): " + "; ".join(out))
np.savez_compressed(f'/Users/adriansager/playground/riemann/.scratch/explore54/line_{kind}_{int(T1)}.npz', t=t.astype(np.float64)[::1], eta=eta.astype(np.float32), logF=np.log(np.abs(F)).astype(np.float32))
print(f"  total {time.time()-t0c:.0f}s")
