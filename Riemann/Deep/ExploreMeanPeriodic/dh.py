# explore14: Davenport-Heilbronn adversary in the FRS template Z_f(s) = s(1-s) Λ(2s)Λ(2s-1)/Λ_f(s)
import mpmath as mp, time
mp.mp.dps = 40
pi = mp.pi; t0 = time.time()
r5 = mp.sqrt(5)
kap = (mp.sqrt(10-2*r5)-2)/(r5-1)
A = [1, kap, -kap, -1]
def f(s):  return 5**(-s)*sum(A[a-1]*mp.zeta(s, mp.mpf(a)/5) for a in range(1,5))
def fp(s): return -mp.log(5)*f(s) + 5**(-s)*sum(A[a-1]*mp.zeta(s, mp.mpf(a)/5, 1) for a in range(1,5))
def Lf(s): return (5/pi)**(s/2)*mp.gamma((s+1)/2)*f(s)
def Lam(s): return pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
s = mp.mpc('0.3','2.0')
print('kappa =', mp.nstr(kap, 15), ' FE check Lf(s)/Lf(1-s) =', mp.nstr(Lf(s)/Lf(1-s), 12))
# real values on (0,1)
vals = [(mp.mpf(k)/40, mp.re(Lf(mp.mpf(k)/40))) for k in range(1,40)]
sc = [v for v in vals if v[1] <= 0]
print('Lf on (0,1): min', mp.nstr(min(v[1] for v in vals),8), 'max', mp.nstr(max(v[1] for v in vals),8), ' nonpositive samples:', len(sc))
L12 = Lf(mp.mpf(1)/2)
a2 = -1/(16*L12)
print('Lf(1/2) =', mp.nstr(L12,12), '  a_{-2}^DH = -1/(16 Lf(1/2)) =', mp.nstr(a2,12))
rho = mp.findroot(f, mp.mpc('0.808517','85.699348'))
print('off-line zero rho* =', mp.nstr(rho, 15), ' |f(rho*)| =', mp.nstr(abs(f(rho)),3))
dLf = (5/pi)**(rho/2)*mp.gamma((rho+1)/2)*fp(rho)
c = rho*(1-rho)*Lam(2*rho)*Lam(2*rho-1)/dLf
print('residue c* =', mp.nstr(c, 8), ' |c*| =', mp.nstr(abs(c), 8), ' log10|c*| =', mp.nstr(mp.log10(abs(c)), 8))
beta = mp.re(rho)
# sign flip of x^{1/2} h(x) = a2*t + 2|c*| e^{(beta-1/2)t} cos(...) + tiny: first t where 2|c*| e^{(beta-1/2)t} = |a2| t
g = lambda t: (beta-mp.mpf(1)/2)*t - mp.log(t) - mp.log(abs(a2)/(2*abs(c)))
tf = mp.findroot(g, 700)
print('t_flip =', mp.nstr(tf, 10), '  x_flip = e^{-t} = 10^', mp.nstr(-tf/mp.log(10), 8))
# same estimate for a hypothetical zeta zero beta+iT: |c| ~ |rho(1-rho)Λ(2ρ)Λ(2ρ-1)/(π^{-ρ/2}Γ(ρ/2)ζ'(ρ))| with |ζ'|~1
for (b,T) in [(0.6, 100), (0.6, 1000), (0.87, 3e12), (0.51, 3e12)]:
    rr = mp.mpc(b, T)
    cc = Lam(2*rr)*Lam(2*rr-1)/(pi**(-rr/2)*mp.gamma(rr/2))   # zeta'(rho) ~ 1
    a2z = mp.mpf('0.0628619872093984')
    gg = lambda t: (b-0.5)*t - mp.log(t) - mp.log(a2z/(2*abs(cc)))
    tt = mp.findroot(gg, max(10, float(-mp.log(abs(cc))/(b-0.5))))
    print('zeta-type zero beta=%s T=%s: log10|c|=%s  t_flip=%s  (x=10^%s)' % (b, T, mp.nstr(mp.log10(abs(cc)),8), mp.nstr(tt,8), mp.nstr(-tt/mp.log(10),8)))
print('time %.1f s' % (time.time()-t0))
