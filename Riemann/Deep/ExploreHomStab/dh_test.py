# explore9 DH test: the Davenport-Heilbronn function f = ((1-ik)/2) L(s,psi) + ((1+ik)/2) L(s,psibar), psi mod 5, psi(2)=i,
# has an off-line zero (0.808517182457 + 85.699348485378 i) but its mean square follows the CFKRS-type recipe:
#   int_0^T |f(1/2+it)|^2 dt ~ ((1+k^2)/2)(4/5) T (log(5T/2pi) + 2 gamma - 1 + log5/2) + 2 Re(c1 conj(c2)(1+eps^2)) L(1,(./5)) T
# Usage: python3 dh_test.py T h
import sys
import mpmath as mp
mp.mp.dps = 15
T = float(sys.argv[1]); h = float(sys.argv[2])
psi = {1: 1, 2: 1j, 4: -1, 3: -1j}
kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
c1 = (1 - 1j * kap) / 2; c2 = (1 + 1j * kap) / 2
def Ls(s):
    hz = {a: mp.zeta(s, mp.mpf(a) / 5) for a in range(1, 5)}
    Lp = sum(psi[a] * hz[a] for a in hz) * mp.power(5, -s)
    Lm = sum(mp.conj(psi[a]) * hz[a] for a in hz) * mp.power(5, -s)
    return Lp, Lm
s0 = mp.mpc(0.808517182457, 85.699348485378)
Lp, Lm = Ls(s0)
print("|f(s0)| at the DH zero:", mp.nstr(abs(c1 * Lp + c2 * Lm), 5), " |L(s0,psi)|:", mp.nstr(abs(Lp), 5))
tau = sum(psi[a] * mp.exp(2j * mp.pi * a / 5) for a in psi)
eps = tau / (1j * mp.sqrt(5))
L1 = 2 * mp.log((1 + mp.sqrt(5)) / 2) / mp.sqrt(5)
cross = 2 * mp.re(c1 * mp.conj(c2) * (1 + eps ** 2)) * L1
print("eps(psi) =", mp.nstr(eps, 10), " cross density =", mp.nstr(cross, 10))
n = int(T / h)
If = 0.0; Ip = 0.0
marks = [100, 200, 400, 800, 1200, 1600, 2000]
for j in range(n + 1):
    t = j * h
    Lp, Lm = Ls(mp.mpc(0.5, t))
    wgt = 0.5 if j in (0, n) else 1.0
    If += wgt * float(abs(c1 * Lp + c2 * Lm) ** 2) * h
    Ip += wgt * float(abs(Lp) ** 2) * h
    if any(abs(t - m) < h / 2 for m in marks):
        lg = float(mp.log(5 * t / (2 * mp.pi)) + 2 * mp.euler - 1 + mp.log(5) / 2)
        predp = 0.8 * t * lg
        predf = float((1 + kap ** 2) / 2) * predp + float(cross) * t
        print(f"T={t:7.1f}  int|f|^2={If:10.3f} recipe={predf:10.3f} rel={If/predf-1:+.4f} | int|L(psi)|^2={Ip:10.3f} formula={predp:10.3f} rel={Ip/predp-1:+.4f}", flush=True)
