# obstruct.py (weildefl): size estimates behind DERIVATION.md sec. 7 (NOT a proof; mpmath, 30 digits).
# For the normalised ground rung v = Phi|_W/||Phi|_W||: boundary value v(L/2); jumps w_n v(L/2) of the coupling beta = Q(v, .) at +-(L/2 - log n);
# high-frequency L^2 mass of beta above angular frequency Omega ~ (2/(pi Omega)) v(L/2)^2 sum_{n<c} w_n^2 (two jumps per n, sin^2 averaged);
# its share in the Schur term ~ that / q, q = log(Omega/(2 pi)) + 2 sigma (generous upper bound of the symbol); Omega needed so that the
# share falls below lambda; minorant drop R(v, v) ~ (2 v(L/2)^2/pi) (log(r/(2 pi)) - 1)/r at r = r_{K2} (all odd k > K2 not in K dropped)
# and its part beyond kscan, ~ (2 v(L/2)^2/pi)(log(r_ks/(2 pi)) - 1)/r_ks, which is known only to O(1) relative accuracy.
import mpmath as mp
mp.mp.dps = 30
def Phi(u):
    u = abs(mp.mpf(u)); z = mp.e**(2*u)
    return mp.fsum((2*mp.pi**2*n**4*mp.e**(9*u/2) - 3*mp.pi*n**2*mp.e**(5*u/2))*mp.e**(-mp.pi*n*n*z) for n in range(1, 12))
def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else 0
    return 0
# lambda: U(c) = Rayleigh quotient upper bound (WeilSimple) for c = 5, 10, 16; certified eps (lower bound) for 20, 25
lam = {5: mp.mpf('1.2218e-17'), 10: mp.mpf('1.9445e-43'), 16: mp.mpf('3.7430e-75'), 20: mp.mpf(2)**-327, 25: mp.mpf(2)**-416}
K2 = {5: 201, 10: 201, 16: 427, 20: 667, 25: 1043}; kscan = {5: 201, 10: 14797, 16: 313858, 20: 7319040, 25: 29076033}
print('c  v(L/2)     sum w_n^2  sigma   beta_hi^2(Om=1e3)  Schur share  lambda     Omega_needed  R(v,v)~     R beyond kscan~')
for c in (5, 10, 16, 20, 25):
    L = mp.log(c); a = L/2
    n2 = 2*mp.quad(lambda u: Phi(u)**2, [0, a/2, a])
    vb = Phi(a)/mp.sqrt(n2)
    sw2 = mp.fsum(mangoldt(n)**2/n for n in range(2, c)); sig = mp.fsum(mangoldt(n)/mp.sqrt(n) for n in range(2, c))
    Om = mp.mpf(1000); bh = 2/(mp.pi*Om)*vb**2*sw2; q = mp.log(Om/(2*mp.pi)) + 2*sig; share = bh/q
    # Omega needed: (2/(pi Om)) vb^2 sw2 / (log(Om/2pi) + 2 sig) = lambda  -> solve
    f = lambda lo: mp.log(2/mp.pi*vb**2*sw2) - lo - mp.log(lo - mp.log(2*mp.pi) + 2*sig) - mp.log(lam[c])   # decreasing in lo
    lo = mp.findroot(f, (mp.mpf(3), mp.mpf(2000)), solver='bisect')
    rK = mp.pi*K2[c]/L; rks = mp.pi*kscan[c]/L
    R = 2*vb**2/mp.pi*(mp.log(rK/(2*mp.pi)) - 1)/rK; Rks = 2*vb**2/mp.pi*(mp.log(rks/(2*mp.pi)) - 1)/rks
    print('%-2d %.3e  %.3f      %.3f   %.2e           %.2e     %.2e   %s       %.1e     %.1e' % (c, vb, sw2, sig, bh, share, lam[c], mp.nstr(mp.e**lo, 2), R, Rks))
