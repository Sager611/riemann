
# Question 5 of 2508.11108 (P=x, beta=1, theta=sqrt(3/5)/R): greedy chain of nodes R_a with zero-extension.
# Condition for all R in [R_a, R_b]:  (2/sqrt15)(I_a + e^{-R_a}) + (1/2 - sqrt(3/5)/3 - 2/sqrt15) e^{-R_b} < e^{-(2/3)sqrt(3/5)}
import numpy as np, sys, json
sys.argv = [sys.argv[0]]  # opt2 reads argv
exec(open('opt2.py').read().split("from fractions import Fraction as Fr")[0])
s15 = np.sqrt(15); a0 = 2/s15; b0 = 0.5 - np.sqrt(0.6)/3 - 2/s15; rhs = np.exp(-(2/3)*np.sqrt(0.6))
def I_of(R, M):
    xi = np.sqrt(0.6)/R
    G0, G1, end = basis(R, M)
    h, c, v = solveQ(G0, G1, end, 1.0/xi, xi/3.0, True)
    # c = 1/2 - sqrt(3/5)/3 + (2/sqrt15) e^R I  =>  I = (c - 1/2 + sqrt(3/5)/3) / ((2/sqrt15) e^R)
    I = (c - 0.5 + np.sqrt(0.6)/3)/(a0*np.exp(R))
    return I, h, c
M = 30
R = 1.5491; nodes = []
while True:
    I, h, c = I_of(R, M)
    lhs = a0*(I + np.exp(-R))
    gap = rhs - lhs   # if > 0: all R' >= R covered
    if gap > 0:
        nodes.append((R, I, None, h)); print(f"R_a={R:.6f} I={I:.8f} covers [R_a, inf), margin {gap:.3e}"); break
    # need -b0 e^{-R_b} > lhs - rhs  i.e. e^{-R_b} > (lhs-rhs)/(-b0)
    eb = (lhs - rhs)/(-b0)
    if eb >= np.exp(-R):
        print(f"FAIL at R_a={R:.6f}: I={I:.8f}, kappa-2theta/3 condition fails even at R_a"); break
    Rb = -np.log(eb)
    Rb_use = R + 0.999*(Rb - R)
    Rb_use = np.floor((R + 0.99*(Rb - R))*1e4)/1e4
    nodes.append((R, I, Rb_use, h))
    R = Rb_use
    if len(nodes) > 20000: print('too many'); break
print(len(nodes), "nodes:", [(a, round(b,6), c) for a,b,c,_ in nodes])
json.dump([{'R': f'{a:.4f}', 'Rb': (None if c is None else f'{c:.4f}'), 'I_float': b, 'h': [repr(float(x)) for x in hh]} for a,b,c,hh in nodes], open('q5nodes.json','w'), indent=1)
