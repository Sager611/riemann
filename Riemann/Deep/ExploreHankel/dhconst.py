# explore13: Davenport-Heilbronn constants for the Nyman-Beurling/Baez-Duarte comparison.
#  f(s) = 5^{-s} sum_r a_r zeta(s, r/5), a = (1, k, -k, -1) for r = 1..4, k = (sqrt(10-2sqrt5)-2)/(sqrt5-1).
#  xi_f(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s) = xi_f(1-s), real on the critical line.
#  C_f = sum_rho 1/(rho(1-rho)) = 2 xi_f'/xi_f(1)   (zeta analogue: 2 + gamma - log 4 pi).
#  Off-line zeros rho (Re > 1/2) give the Beurling floor: lim d_N^2 >= || P_-( 1/(s B(s)) ) ||^2, B = Blaschke product.
import mpmath as mp
mp.mp.dps = 30
kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = {1: mp.mpf(1), 2: kap, 3: -kap, 4: mp.mpf(-1)}

def f(s):
    return mp.power(5, -s) * sum(A[r] * mp.zeta(s, mp.mpf(r) / 5) for r in A)

def xi(s):
    return mp.power(5 / mp.pi, (s + 1) / 2) * mp.gamma((s + 1) / 2) * f(s)

f1 = -sum(A[r] * mp.digamma(mp.mpf(r) / 5) for r in A) / 5
fp1 = mp.diff(f, 1)
C_f = 2 * (mp.log(5 / mp.pi) / 2 - mp.euler / 2 + fp1 / f1)
C_z = 2 + mp.euler - mp.log(4 * mp.pi)
print("kappa =", mp.nstr(kap, 15), " f(1) =", mp.nstr(f1, 15), "(check", mp.nstr(f(mp.mpf(1) + mp.mpf(10)**-12), 12), ")")
print("C_f = sum 1/(rho(1-rho)) =", mp.nstr(C_f, 12), "   C_zeta =", mp.nstr(C_z, 12))

print("first on-line zeros of f (t):")
prev = None; zs = []
t = mp.mpf('0.05')
while t < 40 and len(zs) < 12:
    v = mp.re(xi(mp.mpc(0.5, t)))
    if prev is not None and v * prev < 0:
        z = mp.findroot(lambda u: mp.re(xi(mp.mpc(0.5, u))), (t - mp.mpf('0.05'), t), solver='anderson')
        zs.append(z)
    prev = v; t += mp.mpf('0.05')
print("  ", [mp.nstr(z, 10) for z in zs])
print("  1/|rho_1| =", mp.nstr(1 / abs(mp.mpc(0.5, zs[0])), 8))

guesses = [mp.mpc('0.808517', '85.699348'), mp.mpc('0.650830', '114.163343'),
           mp.mpc('0.574356', '166.479306'), mp.mpc('0.724258', '176.702461')]
off = []
for g in guesses:
    try:
        z = mp.findroot(f, g)
        ok = abs(f(z)) < mp.mpf(10)**-15 and mp.re(z) > 0.5 and abs(z - g) < 0.01
        print("  off-line zero", mp.nstr(z, 14), "|f| =", mp.nstr(abs(f(z)), 3), "accepted" if ok else "REJECTED")
        if ok: off.append(z)
    except Exception as e:
        print("  no convergence near", g, e)

rh = off + [mp.conj(z) for z in off]
c = []
for r in rh:
    rest = mp.mpf(1)
    for q in rh:
        if q != r:
            rest *= (r - q) / (r - (1 - mp.conj(q)))
    c.append((r - (1 - mp.conj(r))) / (r * rest))
floor = mp.re(sum(c[j] * mp.conj(c[l]) / (rh[j] + mp.conj(rh[l]) - 1) for j in range(len(rh)) for l in range(len(rh))))
approx = sum((2 * mp.re(r) - 1) / abs(r)**2 for r in rh)
print("Beurling floor from", len(rh), "off-line zeros: lim d_N^2 >=", mp.nstr(floor, 8), " (sum (2 sigma-1)/|rho|^2 =", mp.nstr(approx, 8), ")")
print("first pair only:", mp.nstr(sum((2 * mp.re(r) - 1) / abs(r)**2 for r in [off[0], mp.conj(off[0])]), 8))
print("log N where C_f/log N = floor:", mp.nstr(C_f / floor, 6), "  i.e. N ~ 10^", mp.nstr(C_f / floor / mp.log(10), 5))
offC = sum(4 * mp.re(1 / (r * (1 - r))) for r in off)
print("off-line part of C_f from these quadruples:", mp.nstr(offC, 6), " -> on-line (Burnol) part >=", mp.nstr(C_f - offC, 10), "(minus unknown higher off-line zeros)")
