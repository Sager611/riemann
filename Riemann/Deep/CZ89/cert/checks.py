
# Arb checks of the scalar claims used by the paper (python-flint, 96-bit balls).
from flint import arb, ctx
import table as T
ctx.prec = 96
PI = arb.pi(); b = arb(1)/arb(2).sqrt()
def fstar(u):   # u may be a ball
    return (arb(2).sqrt()*u).cos()/(arb(2).sqrt()*b.sin()) + arb(2)/25*(2*PI*u).cos() + arb(1)/25*(4*PI*u).cos() - arb(1)/200*(6*PI*u).cos()
# 1) f_* >= 0.79 on [0,1/2] (f_* even): 4096 sub-balls
N = 4096; worst = None
for i in range(N):
    X = arb((2*i+1), 0) / (4*N) ; X = arb(X.mid(), arb(1)/(4*N))
    v = fstar(X)
    lo = v.lower()
    worst = lo if worst is None or lo < worst else worst
print("min_{[0,1/2]} f_* >=", worst.str(10), " > 0.79:", bool(worst > arb('0.79')))
# 2) K_*(0.743), K_*(1/2)
k743 = T.Kders(arb('0.743'))[0]; k12 = T.Kders(arb(1)/2)[0]
CMT = arb(1)/2 + b.cot()/arb(2).sqrt()
k12c = 2*PI*(CMT - arb(1)/2)/(PI*PI - 2) + arb(167)/(3500*PI)
print("K_*(0.743) =", k743.str(15), " > 0.38029:", bool(k743 > arb('0.38029')))
print("K_*(1/2) =", k12.str(15), " closed form", k12c.str(15), " > 2/3:", bool(k12 > arb(2)/3))
print("0.306*K(0.743)^2 - 0.044 =", (arb('0.306')*k743*k743 - arb('0.044')).str(6))
# 3) final constant
C = CMT + arb(321)/80000 - arb(2449)/(1440000*PI*PI)
tau = 2 + arb(4791).sqrt()/50
U = 1 - 4*(C - arb(10133)/10000)/(tau*tau)
print("C_MT =", CMT.str(15), " C =", C.str(15))
print("U =", U.str(15), " U - 8/9 =", (U - arb(8)/9).str(6), " >0:", bool(U > arb(8)/9))
alpha = arb(209)/5000
print("tau = 2+sqrt(2-2alpha):", bool(abs(tau - 2 - (2-2*alpha).sqrt()) < arb('1e-25')), " (tau-1)^2 =", ((tau-1)**2).str(10), " 114*11/250 =", (arb(114)*11/250).str(6))
# 4) S at a rational point near the numerical minimiser: upper bound for min S
rows = {0:['0.306','0.645','1','1','1','2'],1:['0.347','0.047','0','0','1'],2:['0.347','0.616','0','1'],3:['0.347','0.047','1'],4:['0.347','0.645'],5:['0.306']}
W = {}
for i,v in rows.items():
    for t,val in enumerate(v): W[(i,i+1+t)] = arb(val)
def S(g):
    y = [arb(0)]
    for gi in g: y.append(y[-1] + gi)
    tot = sum((W[(i,j)]*T.Kders(y[j]-y[i])[0]**2 for (i,j) in W), arb(0))
    return tot + sum(g, arb(0))/200
gstar = [arb(x) for x in ['1.0837079364','1.8939190725','1.0413580806','1.0248802119','1.0289099949','1.0255587059']]
Sg = S(gstar)
print("S(g*) =", Sg.str(15), " S(g*) - 11/250 =", (Sg - arb(11)/250).str(6))
g2 = [arb(x) for x in ['1.0426172565','1.0469177998','1.9030981233','1.0568326713','1.0297419287','1.0273947729']]
print("S(g2) =", S(g2).str(15))

