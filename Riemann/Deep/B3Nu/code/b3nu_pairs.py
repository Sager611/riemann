# b3nu: exponent pairs (k,l) maximizing the gain of Lemma C:  nu1(theta) = [theta(1+2l) - (k+l)] / (2(1+l-k));
# gain over (3theta-1)/4 equals (1-theta)(1-l-3k)/(4(1+l-k)).  Pairs generated from (0,1), Bourgain 2017 (13/84,55/84)+eps,
# Huxley 2005 (32/205,269/410)+eps by the A and B processes (A(k,l) = (k/(2k+2),(k+l+1)/(2k+2)), B(k,l) = (l-1/2,k+1/2)).
from fractions import Fraction as F
A = lambda p: (p[0] / (2 * p[0] + 2), (p[0] + p[1] + 1) / (2 * p[0] + 2))
B = lambda p: (p[1] - F(1, 2), p[0] + F(1, 2))
seeds = {'(0,1)': (F(0), F(1)), 'Bourgain(13/84,55/84)+eps': (F(13, 84), F(55, 84)), 'Huxley(32/205,269/410)+eps': (F(32, 205), F(269, 410))}
pairs = {}
for name, s in seeds.items():
    front = [(name, s)]
    for depth in range(9):
        new = []
        for (nm, p) in front:
            for lab, op in (('A', A), ('B', B)):
                q = op(p)
                if 0 <= q[0] <= F(1, 2) <= q[1] <= 1:
                    key = q
                    if key not in pairs: pairs[key] = lab + nm; new.append((lab + nm, q))
        front = new
obj = lambda p: (1 - p[1] - 3 * p[0]) / (1 + p[1] - p[0])
best = sorted(pairs.items(), key=lambda kv: -obj(kv[0]))[:8]
th = F(4947, 10000)
nu1 = lambda p, th: (th * (1 + 2 * p[1]) - (p[0] + p[1])) / (2 * (1 + p[1] - p[0]))
print("top pairs by (1-l-3k)/(1+l-k):")
for p, nm in best:
    print(f"  {nm:40s} (k,l)=({p[0]},{p[1]}) = ({float(p[0]):.6f},{float(p[1]):.6f})  obj={float(obj(p)):.6f}  nu1(0.4947)={float(nu1(p, th)):.6f}  (3th-1)/4={float((3*th-1)/4):.6f}")
p = (F(1, 30), F(13, 15)); print(f"classical A^3B(0,1)=(1/30,13/15): obj={float(obj(p)):.6f} nu1(th)=(82th-27)/110 check: {nu1(p, th) == (82*th-27)/110}  nu1(1)={nu1(p, F(1))}  th0(nu1=0)={(p[0]+p[1])/(1+2*p[1])}")

# classical tree only (seed (0,1)), outermost process A (linear-phase immunity)
cl = {}
front = [('', (F(0), F(1)))]
for depth in range(12):
    new = []
    for (nm, p) in front:
        for lab, op in (('A', A), ('B', B)):
            q = op(p)
            if 0 <= q[0] <= F(1, 2) <= q[1] <= 1 and q not in cl:
                cl[q] = lab + nm; new.append((lab + nm, q))
    front = new
bestc = sorted([(p, nm) for p, nm in cl.items() if nm.startswith('A')], key=lambda t: -obj(t[0]))[:4]
for p, nm in bestc:
    print(f"classical {nm+'(0,1)':16s} (k,l)=({p[0]},{p[1]}) obj={float(obj(p)):.6f} gain/(1-theta) = {float(obj(p))/4:.6f}")
