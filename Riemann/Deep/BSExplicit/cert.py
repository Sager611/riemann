# bsexpl certificate: every constant of DERIVATION.md in arb balls (python-flint), 200-bit precision.
# Run: OMP_NUM_THREADS=1 uv run --no-project --with python-flint python3 cert.py > cert_output.txt   (< 2 s, < 50 MB)
import sys
from flint import arb, ctx

ctx.prec = 200
ok = True


def claim(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


def up(x, d=6):   # printable upper end
    return f"{float((x.mid() + x.rad())):.{d}g}"


def lo(x, d=6):   # printable lower end
    return f"{float((x.mid() - x.rad())):.{d}g}"


s2 = arb(2).sqrt()
D = 6 - 4 * s2                        # Lemma D constant
cq = D / 4                            # quadratic k: (3 - 2 sqrt 2)/2
B2 = (777 - arb(1081).sqrt()) / 852   # Liu's contact B_2
Tstar = B2 + arb("2.2e-254")          # Lean-checked half-plane (Bebel, 1 Nov 2026), for every Dirichlet L-function
log3, log9 = arb(3).log(), arb(9).log()
print("== §1 constants")
print("6 - 4 sqrt2 in [", lo(D, 12), ",", up(D, 12), "]")
print("(3 - 2 sqrt2)/2 in [", lo(cq, 12), ",", up(cq, 12), "]")
print("B_2 =", lo(B2, 15), "; 1 - T* >=", lo(1 - Tstar, 12))
hp9 = (1 - Tstar) * log9
hp3 = (1 - Tstar) * log3
print("(1 - T*) log 9 >=", lo(hp9, 8), "; (1 - T*) log 3 >=", lo(hp3, 8))
claim("(1-T*) log 9 > (3-2sqrt2)/2  [case (a) never binds for n >= 2]", hp9 > cq)
claim("(1-7/8) log 9 > (3-2sqrt2)/2 [same with the 7/8 Lean port]", (arb(1) / 8) * log9 > cq)

print("== §2 Lemma D: b <= 6-4sqrt2 => contradiction at a = (1+sqrt2) b")
for bnum in ["0.001", "0.1", "0.3", "0.343"]:
    b = arb(bnum)
    a = (1 + s2) * b
    lhs = 2 / (a + b) - 1 / a
    claim(f"b={bnum}: 2/(a+b)-1/a = (3-2sqrt2)/b", abs(lhs - (3 - 2 * s2) / b) < arb("1e-50"))
    claim(f"b={bnum}: (3-2sqrt2)/b >= 1/2 (contradiction with 2/(s-beta) < 1/(s-1) + L/2)", (3 - 2 * s2) / b >= arb(1) / 2)
amax = (1 + s2) * D                   # = 2 sqrt2 - 2
claim("a_max = (1+sqrt2)(6-4sqrt2) = 2sqrt2-2 < log 3 <= log d_M  (sigma < 2)", (amax < log3) and abs(amax - (2 * s2 - 2)) < arb("1e-50"))

print("== §3 Stark Lemma 4: c4 = 2 exp(21/8 + (gamma + 2 log 2)/(4 log 3))")
psi_half = arb("0.5").digamma()       # = -gamma - 2 log 2
c2 = 2 / log3
c3 = (arb(21) / 8 - (c2 / 8) * psi_half).exp()
c4 = 2 * c3
print("psi(1/2) =", lo(psi_half, 12), "; c3 in [", lo(c3, 10), ",", up(c3, 10), "]; c4 <=", up(c4, 10))
claim("c4 < 43.166", c4 < arb("43.166"))
# c2 = 2/log 3 needs d_K >= 3^{n/2}: n = 2 trivial (d >= 3); n >= 3 by Minkowski d >= (pi/4)^n (n^n/n!)^2 and
# n! <= e sqrt(n) (n/e)^n, i.e. d >= (pi e^2/4)^n/(e^2 n); need (pi e^2/(4 sqrt3))^n >= e^2 n for n >= 3.
pi, e = arb.pi(), arb(1).exp()
q = pi * e * e / (4 * arb(3).sqrt())
claim("q = pi e^2/(4 sqrt3) > 3.35", q > arb("3.35"))
claim("q^3 >= 3 e^2 (n = 3), and q > (n+1)/n for n >= 1 so q^n/n increases", (q ** 3 >= 3 * e * e) and (q > 2))

print("== §4 Theorem A: c(k) = (6 - 4 sqrt2)/(2^{n-1} n!), n = [k:Q]; quadratic k: (3 - 2sqrt2)/2")
fact = 1
for n in range(2, 9):
    fact *= n
    lam = 2 ** (n - 1) * fact
    c = D / lam
    guo = min(n * 0.01, float(lo(c, 12)))
    print(f"n={n}: lambda <= {lam:>8d}, c(k) >= {lo(c, 6)}   (with Guo's c_Q = 1/100 instead of the half-plane: >= {guo:.6g})")
claim("quadratic k: c(k) = (3-2sqrt2)/2 >= 0.085786", cq >= arb("0.085786"))

print("== §5 Theorem B: h_K R_K > c_BS d_K^{1/2} / ((2 pi)^{[K:Q]/2} log d_K), c_BS = 2 min(1/4, c_1)/c4")
fact = 1
for n in range(2, 7):
    fact *= n
    c = cq if n == 2 else D / (2 ** (n - 1) * fact)
    eta = min(float(lo(hp3, 12)), float(lo(D / fact, 12)))
    c1 = min(float(lo(D, 12)), float(lo(c, 12)), eta)
    cbs = 2 * min(0.25, c1) / float(up(c4, 12))
    print(f"n={n}: c(k) >= {float(lo(c, 12)):.6g}, eta_k >= {eta:.6g}, c_1 >= {c1:.6g}, c_BS >= {cbs:.6g}")
cbs2 = 2 * cq / c4
print("quadratic k (Q(i), Q(sqrt-3), ...): c_BS >=", lo(cbs2, 8))
claim("quadratic k: c_BS >= 3.97e-3", cbs2 > arb("3.97e-3"))
cbsQ = 2 * hp3 / c4
print("k = Q (Galois K/Q, rem:almost): c_1 = min(6-4sqrt2, (1-T*) log 3) >=", lo(hp3, 8), "; c_BS >=", lo(cbsQ, 8))

print("== §6 comparison: general-field half-planes (not used in Thms A, B)")
for name, w in [("Cureton-Thorner 31/32 (Lean)", arb(1) / 32), ("Liu 2610.26835 B=0.874886647040", 1 - arb("0.874886647040"))]:
    print(f"{name}: (1-beta) log d_F >= {lo(w, 6)} log d_F; Q(i) (d_F >= 144): {lo(w * arb(144).log(), 6)}; Q(sqrt-3) (d_F >= 117): {lo(w * arb(117).log(), 6)}")
cbsCT = 2 * (log3 / 32) / c4
print("CT for every K != Q: c_BS = 2 min(1/4, log3/32)/c4 >=", lo(cbsCT, 6))

print("== §7 k = Q specialization: Theorem A with M = F gives min(c_Q, 6-4sqrt2) <= c_Q (no gain over the input)")
for name, cQ in [("Guo 2610.30476", arb("0.01")), ("2610.21511", arb("5.811e-4")), ("Lean half-plane (1-T*) log 3", hp3)]:
    claim(f"{name}: min(c_Q, 6-4sqrt2) = c_Q", (cQ < D))
print("ALL PASS" if ok else "SOME FAIL")
sys.exit(0 if ok else 1)
