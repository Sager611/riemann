# pyright: reportMissingImports=false
"""Rigorous Arb enclosures of gamma_k = int_0^oo Phi(u) u^{2k}/(2k)! du (Gershon's moment sequence, any
u-scaling gives the same L_r, Theta_r), and of C_1(n) = n^2 |log Theta_1(n)| (Gershon, Lemma 11.2/11.6):
  L1(n) = g_n^2/(g_{n-1} g_{n+1}),  eps(n) = 1 - 1/L1(n),  Theta1(n) = eps(n)^2/(eps(n-1) eps(n+1)).
Phi(u) = sum_{m>=1} (2 pi^2 m^4 e^{9u} - 3 pi m^2 e^{5u}) exp(-pi m^2 e^{4u}).
Integral over [0,U] by acb.integral with terms m<=M; rigorous tail bound added as ball radius:
 with t = e^{4u} >= 1, u <= t/4, e^{9u} = t^{9/4}, du = dt/(4t):
 |term_m(u)| u^{2k} du <= 2 pi^2 m^4 4^{-2k-1} t^{2k+5/4} e^{-pi m^2 t} dt  (also bounds the 3pi m^2 e^{5u} part:
 2pi^2 m^4 e^{9u} >= 3 pi m^2 e^{5u}, so |term| <= 2pi^2 m^4 e^{9u} e^{-pi m^2 t}).
 m > M, all u:  sum <= sum_m 2 pi^2 m^4 4^{-2k-1} Gamma(a+1, pi m^2)/(pi m^2)^{a+1}, a = 2k+5/4 (summed to m=M+2000, + geometric-ish tail bound).
 m <= M, u > U: 2 pi^2 m^4 4^{-2k-1} Gamma(a+1, pi m^2 T)/(pi m^2)^{a+1}, T = e^{4U}."""
import sys
from flint import arb, acb, ctx
ctx.prec = 400
M, U = 6, arb(3)/2
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 54
pi = arb.pi()

def tail(k):
    a = arb(2*k) + arb(5)/4
    T = (4*U).exp()
    tot = arb(0)
    for m in range(1, M+1):
        c = pi*m*m
        tot += 2*pi**2*m**4 / arb(4)**(2*k+1) * (c*T).gamma_upper(a+1) / c**(a+1)  # flint: z.gamma_upper(s) = Gamma(s, z)
    # m > M: explicit sum m = M+1..M+200; remainder m >= M+201 = 207 rigorous: for c = pi m^2 >= 2(a+1) and t >= 1,
    # t^a e^{-ct} <= e^{-ct/2} (t^a e^{-ct/2} decreasing on t >= 1), so int_1^oo <= (2/c) e^{-c/2};
    # then 2 pi^2 m^4 4^{-2k-1} (2/c) e^{-c/2} <= 4 pi m^2 e^{-pi m^2/2} <= e^{-m} (m >= 2), sum_{m>=207} e^{-m} < 2 e^{-207}.
    for m in range(M+1, M+201):
        c = pi*m*m
        tot += 2*pi**2*m**4 / arb(4)**(2*k+1) * c.gamma_upper(a+1) / c**(a+1)  # t >= 1
    assert pi*(M+201)**2 >= 2*(a+1)
    tot += 2*arb(-207).exp()
    return tot

def gamma_k(k):
    def f(u, _):
        s = acb(0)
        e4 = (4*u).exp()
        for m in range(1, M+1):
            s += (2*pi**2*m**4*(9*u).exp() - 3*pi*m**2*(5*u).exp()) * (-pi*m*m*e4).exp()
        return s * u**(2*k)
    I = acb.integral(f, 0, U).real
    t = tail(k)
    I = I + arb(0, t.upper())
    return I / arb(2*k+1).gamma()

g = [gamma_k(k) for k in range(KMAX+1)]
L1 = {n: g[n]**2/(g[n-1]*g[n+1]) for n in range(1, KMAX)}
eps = {n: 1 - 1/L1[n] for n in L1}
C1 = {}
for n in range(2, KMAX-1):
    Th = eps[n]**2/(eps[n-1]*eps[n+1])
    C1[n] = -(n*n)*Th.log()
for n in sorted(C1):
    if n in (10, 20, 30, 40, 50) or n == max(C1):
        print(n, "L1-1 =", (L1[n]-1).str(12, radius=False), " C1 =", C1[n].str(20), " Theta1<1:", bool(C1[n] > 0))
if 50 in C1:
    d = C1[50] - C1[30]
    print("C1(50) - C1(30) in", d.str(20))
    print("certified: C1(50)-C1(30) > 3/100 :", bool(d > arb(3)/100))
    print("Lemma 11.6 (sum K_s < 2, delta = 0.452) allows at most", (2*(arb(452)/1000)**30 + 2*(arb(452)/1000)**50).str(5))
print("rel radius gamma_KMAX:", (g[KMAX].rad()/g[KMAX].mid()))
