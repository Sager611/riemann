
import json, math
tab = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/BellottiFix/code/tables.json'))
def lT(s):
    s = s.strip(); return math.log(3e12) if s.startswith('3') else float(s[5:-1])
def nulog(L, L0):  # nu(T)*log T with branch chosen by T0 (as in the paper)
    if L0 < 46.2: return 1/5.558691
    if L0 < 170.2:
        J = L/6 + math.log(L) + math.log(0.618)
        return (0.04962 - 0.0196/(J+1.15))/(J + 0.685 + 0.155*math.log(L))*L
    if L0 < 481958: return math.log(L)/21.233
    return 1/(53.989*L**(2/3)*math.log(L)**(1/3))*L
out = []
for i, (cb, p) in enumerate(zip(tab['CB'], tab['uxvw'])):
    L0, L1 = lT(cb[0]), lT(cb[1]); x = float(p[3])
    a = min(2*x*nulog(L, L0) for L in [L0 + (L1-L0)*k/200 for k in range(201)])
    f = math.exp(a-1)/a if a < 1 else 1.0
    out.append((i+1, x, round(a,3), round(f,4)))
print([o for o in out if o[2] < 1.0])
print('min a', min(o[2] for o in out))
