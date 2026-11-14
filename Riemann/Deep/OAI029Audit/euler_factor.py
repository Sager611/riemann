# 029 Prop 7.3: series (05-poisson eq:exact-Euler-factor, C from eq:C-local-values) vs closed form eq:Euler-good-exact,
# and the E1=0 specialization formula for H-1 (l.665-672). Degree-one prime p = 13, 37; random theta, eta.
import cmath, random, sys
from euler_check import setup, T
PLANT = len(sys.argv) > 1 and sys.argv[1] == 'plant'
random.seed(1)
def run(p, xi, w, z, eta, th):
    chi = setup(p); Q = p
    tau = {v: T(p, chi, v, 1) for v in range(6)}
    def Tv(v, y): return (Q-1 if v % 6 == 0 else 0) if y % p == 0 else (chi(y, v).conjugate()*tau[v % 6])
    def C(t, k, J):
        if t == 0: return 1 if k == 0 else (1 if J == 0 else 0)
        if k == 0: return Q**(t-1)*Tv(t, p**(J-t+1)) if J >= t-1 else 0
        if k == 1: return tau[1]*Q**(t-1)*Tv(t-1, p**(J-t)) if J >= t else 0
        return 0
    g1, g2 = tau[1]/Q**0.5, tau[2]/Q**0.5
    vt = -eta/(g1*g2)                                  # gamma1 gamma2 vartheta = -eta
    if PLANT: vt = eta/(g1*g2)                         # planted sign error
    S = 0
    for e0 in (0, 1):
        for l in range(0, 6):
            for k in range(0, 14):
                for mp in range(0, 8):
                    t = e0 + 3*l; J = 6*mp
                    c = C(t, k, J)
                    if c == 0: continue
                    S += (g2**e0)*(vt**t)*(th**(k-t))*c*Q**(-e0*(0.5+xi)-l*(1+3*xi)-k*w-6*mp*z)
    V = Q**(-6*z); W = th*Q**(-w); D1 = eta*th.conjugate()*Q**(-xi); E1 = vt**6*Q**(4-6*xi)*V
    Fc = 1/(1-V) + W/(1-W) + E1/(1-E1)*((1-1/Q)/(1-V) + W) - D1/(1-E1)*(1 + (Q-1)*W*V/(1-V))
    F0 = 1/(1-V) + W/(1-W) - D1*(1 + (Q-1)*W*V/(1-V))
    H0m1 = (1-V)*(1-W)/(1-D1)*F0 - 1
    H0m1_paper = (D1*(W+V-W*V) - W*V - (Q-1)*D1*(1-W)*W*V)/(1-D1)
    return abs(S-Fc)/abs(Fc), abs(H0m1-H0m1_paper)/max(abs(H0m1), 1e-300)
worst1 = worst2 = 0
for p in (13, 37):
    for trial in range(6):
        xi = complex(1.2+random.random(), 3*random.random()); w = complex(1.1+random.random(), random.random())
        z = complex(0.35+random.random(), random.random()); eta = cmath.exp(2j*cmath.pi*random.random())
        th = cmath.exp(2j*cmath.pi*random.randrange(6)/6)
        a, b = run(p, xi, w, z, eta, th); worst1 = max(worst1, a); worst2 = max(worst2, b)
print('Euler closed form vs series: max rel diff %.2e ; H(E1=0)-1 formula: max rel diff %.2e' % (worst1, worst2))
