# FodgeAudit consts2.py: constants of Lemma C.5 as used in Lemmas C.6 (N >= 5e4, 2nd differences) and C.8 (N >= 2e4, 3rd differences), Arb 256 bits;
# exactness of the parabola-weight closed forms (parab_forms.py of the bundle) used by the Prop. C.7 code.
import sys, random
from flint import arb, ctx
ctx.prec = 256
pi = arb.pi(); z3 = arb(3).zeta(); half = arb(1)/2
def c5(N, p, k3, k4):      # h = ((1+1/75) N/N')^p, eps3 = k3/N', eps4 = k4/N'^2, b - a <= L <= N/75
    Np = arb(N) - half; h = ((1 + arb(1)/75)*arb(N)/Np)**p; e3 = k3/Np; e4 = k4/Np**2; L = arb(N)/75
    X = e3**2*(h/6 + h**1.5/2)**2 + h.sqrt()*(e3**2*(7*h*h/36 + h**2.5/2) + e4*(h/12 + h**1.5/3))
    return h, 8/pi.sqrt() + (h.sqrt()-1)/2, 1 + e3/8 + h**2.5*X*L*L/(4*pi), 7*z3/pi**2*h**1.5 + 21*z3/(2*pi**2)*h
for lab, N, p, k3, k4, claims, used in (('C.6 (2nd diff.)', 5*10**4, 4, 4, 20, ('4.52695', '1.00051', '2.27161'), ('4.527', '1.0006', '2.272')),
                                        ('C.8 (3rd diff.)', 2*10**4, 5, 5, 30, ('4.53039', '1.00085', '2.30835'), ('4.531', '1.001', '2.31')),
                                        ('C.3 (1st diff.)', 10**4, 3, 3, 12, ('4.52359', '1.00031', '2.23600'), ('4.524', '1.0004', '2.24'))):
    h, K, C0, E1 = c5(N, p, k3, k4)
    print(f'{lab}: h <= {h.upper().str(10, radius=False)}')
    for nm, v, cl, us in zip(('8/sqrt(pi)+(sqrt h-1)/2', '1+eps3/8+h^2.5 Xi L^2/(4pi)', '7z3h^1.5/pi^2+21z3h/(2pi^2)'), (K, C0, E1), claims, used):
        print(f'   {nm} = {v.upper().str(12, radius=False)}  <= printed {cl}: {v < arb(cl)} (gap {(arb(cl)-v).upper().str(3, radius=False)}),  <= used {us}: {v < arb(us)}')
    # monotone in N: worst case at the smallest N (check at 2N)
    print(f'   at 2N: {[x.upper().str(8, radius=False) for x in c5(2*N, p, k3, k4)[1:]]}')
sys.path.insert(0, sys.argv[1])
from parab_forms import W1_num, W1_den, r_num, r_den, R_num, R_den
random.seed(3); bad = 0; n = 0
for M in list(range(1, 60)) + [random.randint(60, 400) for _ in range(40)]:
    w = [(j+1)*(M-j) for j in range(M)]
    r = [sum(w[j]*w[j+m] for j in range(M-m)) for m in range(M)]
    if W1_num(M) != W1_den*sum(w): bad += 1
    for m in (list(range(M)) if M < 60 else random.sample(range(M), 20)):
        n += 1
        if r_num(M, m) != r_den*r[m]: bad += 1
    for _ in range(20):
        a = random.randint(1, max(1, M-1)); b = random.randint(a, max(a, M-1))
        if b <= M-1 and R_num(M, b) - R_num(M, a-1) != R_den*sum(r[a:b+1]): bad += 1
print(f'parab_forms: W1, r(m), group sums sum_(m=a..b) r(m) exact for M < 60 (all m) and 40 random M in [60, 400]: {bad} mismatches ({n} r(m) checks)')
