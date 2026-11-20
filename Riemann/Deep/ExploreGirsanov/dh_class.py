# dh_class.py (explore56): adversarial check. Is Davenport-Heilbronn in the positive symmetric class (S2), and does S3 apply?
# DH (as ExploreHelio/dh.py): f(s) = 5^{-s} sum_r a_r zeta(s, r/5), a = (1, k, -k, -1), k = kappa; odd, conductor 5.
# Lambda(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s) = int_0^inf theta(x) x^{(s+1)/2} dx/x, theta(x) = sum_n c_n n exp(-pi n^2 x/5).
# Positive law <=> theta > 0 on (0, inf). Then F(s) = Lambda(s)/Lambda(1) = E[Y^s] with F(1-s) = F(s).
import mpmath as mp
mp.mp.dps = 30
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
a = {1: mp.mpf(1), 2: k, 3: -k, 4: mp.mpf(-1), 0: mp.mpf(0)}
def theta(x, N=400):
    return sum(a[n % 5]*n*mp.e**(-mp.pi*n*n*x/5) for n in range(1, N+1))
xs = [mp.mpf(10)**(mp.mpf(j)/20) for j in range(-40, 41)]   # x in [0.01, 100]
vals = [theta(x) for x in xs]
i = min(range(len(xs)), key=lambda j: vals[j]/abs(theta(xs[j], 1)))
print('kappa =', mp.nstr(k, 15))
print('theta_DH > 0 on the grid 0.01..100:', all(v > 0 for v in vals), ' min theta/first-term =', mp.nstr(vals[i]/theta(xs[i], 1), 6), 'at x =', mp.nstr(xs[i], 4))
print('self-dual check theta(1/x) / (x^{3/2} theta(x)) at x = 2, 3:', [mp.nstr(theta(1/mp.mpf(x))/(mp.mpf(x)**1.5*theta(mp.mpf(x))), 12) for x in (2, 3)])
def Lam(s):
    return (5/mp.pi)**((s+1)/2)*mp.gamma((s+1)/2)*5**(-s)*sum(a[r]*mp.zeta(s, mp.mpf(r)/5) for r in (1, 2, 3, 4))
print('FE check Lambda(0.3+2i) - Lambda(0.7-2i):', mp.nstr(abs(Lam(mp.mpc(0.3, 2)) - Lam(mp.mpc(0.7, -2))), 3))
L0 = Lam(mp.mpf(0))  # = Lambda(1) by the FE; s = 1 itself is a removable 0*inf point of the Hurwitz sum
kt = [mp.diff(lambda s: mp.log(Lam(s)/L0), 1, m, method='quad', radius=0.25) for m in (1, 2, 3)]
l1 = kt[0]; l2 = 2*kt[0] + kt[1]; l3 = 3*kt[0] + 3*kt[1] + kt[2]/2
print('DH Esscher cumulants kt1..3 =', [mp.nstr(x, 10) for x in kt])
print('DH lambda_1, lambda_2, lambda_3 =', mp.nstr(l1, 10), mp.nstr(l2, 10), mp.nstr(l3, 10))
print('S3 lower bound 3k1 - k1^3/2 + 1.5 k2 (2 - k1) =', mp.nstr(3*kt[0] - kt[0]**3/2 + 1.5*kt[1]*(2 - kt[0]), 10))
