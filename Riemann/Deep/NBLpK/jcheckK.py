# nflp: numeric check of the residue identity in A4(b) of PROOF.md for zeta_K, K = Q(i), Q(sqrt 5).
#   J := int_R R(s')/|s'|^2 dt, s' = sigma' + it, R(s) = zeta_K(s)/zeta_K(s+eps), sigma' = 1/p - eps/2;
#   formula: J = 2 pi [ R(2 sigma')/(2 sigma') - c_K / ((2 sigma' - 1) zeta_K(1+eps)) ].
# Also B(eps)/eps, B(eps) = (1/2pi) int |1 - R|^2/|s|^2 dt (the quantity bounded by K_{p,K} eps in A4).
# zeta_K = zeta * L(., chi_D).  Simpson on [0, T] (R(conj s) = conj R(s)), step 0.0005 on [0,2] (pole of zeta_K near t = 0 when
# sigma' is close to 1) and 0.025 on [2,T]; tail 2 int_T^inf dt/(sigma'^2+t^2) (mean Re R = 1).
# Run: perl -e 'alarm 1800; exec @ARGV' uv run --no-project --with mpmath python3 jcheckK.py [field p eps T]   (output: jcheckK.txt)
import mpmath as mp
mp.mp.dps = 15
FIELDS = {"Q(i)": ([0, 1, 0, -1], mp.pi / 4), "Q(sqrt5)": ([0, 1, -1, -1, 1], 2 * mp.log((1 + mp.sqrt(5)) / 2) / mp.sqrt(5))}

def zK(s, chi):
    return mp.zeta(s) * mp.dirichlet(s, chi)

def simpson(f, a, b, h):
    n = int(round((b - a) / h)); n += n % 2; h = (b - a) / n
    S = mp.mpf(0); SB = mp.mpf(0)
    for i in range(n + 1):
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        x, y = f(a + i * h); S += w * x; SB += w * y
    return S * h / 3, SB * h / 3

def run(name, chi, cK, p, eps, T=100, h=mp.mpf("0.025"), tf=2, hf=mp.mpf("0.0005")):
    # fine step on [0, tf]: the pole of zeta_K at s = 1 lies at distance 1 - sigma' (0.046 for p = 32/31, eps = 0.03)
    sp = mp.mpf(1) / p - eps / 2
    def f(t):
        s = mp.mpc(sp, t)
        R = zK(s, chi) / zK(s + eps, chi)
        den = sp * sp + t * t
        return mp.re(R) / den, abs(1 - R) ** 2 / den
    a1, b1 = simpson(f, mp.mpf(0), mp.mpf(tf), hf)
    a2, b2 = simpson(f, mp.mpf(tf), mp.mpf(T), h)
    SJ, SB = a1 + a2, b1 + b2
    tail = (mp.pi / 2 - mp.atan(T / sp)) / sp
    Jnum = 2 * (SJ + tail)
    Bnum = 2 * SB / (2 * mp.pi)  # tail of B not added
    u = 2 * sp
    Jform = 2 * mp.pi * (zK(u, chi) / zK(u + eps, chi) / u - cK / ((u - 1) * zK(1 + eps, chi)))
    Jflip = 2 * mp.pi * (-zK(u, chi) / zK(u + eps, chi) / u - cK / ((u - 1) * zK(1 + eps, chi)))
    Jnopole = 2 * mp.pi * (zK(u, chi) / zK(u + eps, chi) / u)
    print(f"{name} p={mp.nstr(p,6)} eps={mp.nstr(eps,4)} sigma'={mp.nstr(sp,6)}: J_num={mp.nstr(Jnum,9)} J_formula={mp.nstr(Jform,9)}"
          f" diff={mp.nstr(Jnum-Jform,3)} | wrong-sign={mp.nstr(Jflip,7)} no-pole-term={mp.nstr(Jnopole,7)} | B/eps~{mp.nstr(Bnum/eps,5)}"
          f" | c_K/(eps zeta_K(1+eps))={mp.nstr(cK/(eps*zK(1+eps,chi)),6)}", flush=True)

import sys
if len(sys.argv) > 1:  # single case: jcheckK.py <field> <p as a/b> <eps> <T>
    chi, cK = FIELDS[sys.argv[1]]
    a, b = sys.argv[2].split('/')
    run(sys.argv[1], chi, cK, mp.mpf(a) / mp.mpf(b), mp.mpf(sys.argv[3]), T=int(sys.argv[4]))
else:
    for name, (chi, cK) in FIELDS.items():
        run(name, chi, cK, mp.mpf(8) / 7, mp.mpf("0.1"))
        run(name, chi, cK, mp.mpf(32) / 31, mp.mpf("0.03"))
