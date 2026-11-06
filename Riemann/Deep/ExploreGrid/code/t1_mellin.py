# t1_mellin.py (explore21): what hexagonal optimality can say about complex s. With D = theta_L - theta_H >= 0 (Montgomery 1988,
# covolume 1), M(s) = pi^-s Gamma(s) (Z_L(s) - Z_H(s)) = int_0^inf D(a) a^{s-1} da is entire and |M(s)| <= M(sigma).
# At a zero of Z_L or Z_H this is the only constraint; print both sides (L = 2x^2+xy+2y^2, H = hexagonal).
import mpmath as mp
mp.mp.dps = 30
c15 = [0,1,1,0,1,0,0,-1,1,0,0,-1,0,-1,-1]; c3 = [0,1,-1]; c5 = [0,1,-1,-1,1]
def ZL(s): return (mp.sqrt(15)/2)**s*(mp.zeta(s)*mp.dirichlet(s, c15) - mp.dirichlet(s, c3)*mp.dirichlet(s, c5))
def ZH(s): return (mp.sqrt(3)/2)**s*6*mp.zeta(s)*mp.dirichlet(s, c3)
def M(s): return mp.power(mp.pi, -s)*mp.gamma(s)*(ZL(s) - ZH(s))
for sig in (0.5, 0.75807089343824137875, 2):
    print('M(%.5f) = %s   (Z_L - Z_H = %s)' % (sig, mp.nstr(M(sig), 8), mp.nstr(ZL(sig) - ZH(sig), 8)))
s0 = mp.mpc('0.75807089343824137875', '24.482821394256568894')
r1 = mp.findroot(lambda s: mp.dirichlet(s, c3), mp.mpc(0.5, 8.04))
for name, s in (('zero of Z_L (off-line)', s0), ('first zero of L(s,chi_-3) (zero of Z_H)', r1)):
    print('%s s = %s: |M(s)| = %s <= M(sigma) = %s ; |pi^-s Gamma(s)| = %s' % (name, mp.nstr(s, 12), mp.nstr(abs(M(s)), 5), mp.nstr(M(mp.re(s)), 5), mp.nstr(abs(mp.power(mp.pi, -s)*mp.gamma(s)), 5)))
