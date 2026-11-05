# T1b: given positivity, what P(bad) does the PROVABLE first moment give?  P_bad(alpha) <= Zup / h(i alpha)
import json, math
rows=json.load(open('t1_rows.json'))
def h_i(alpha,R): x=alpha*R/2; return R*(math.sinh(x)/x)**2
print("q, A, R, Zup(provable), Z(actual), P_bad<=Zup/h(i a) for a=0.05,0.1,0.25 ; a needed for P<0.1 ; box half-height pi/(3R)")
for r in rows:
    if abs(r['A']-2.0)>1e-9 and not (r['q']==100003 and abs(r['A']-1.8)<1e-9): continue
    R=r['R']; out=[]
    for a in (0.05,0.1,0.25): out.append(r['Zup']/h_i(a,R))
    lo,hi=1e-4,0.5
    for _ in range(60):
        m=(lo+hi)/2
        if r['Zup']/h_i(m,R)<0.1: hi=m
        else: lo=m
    print(f"q={r['q']:6d} A={r['A']:.2f} R={R:6.2f} Zup={r['Zup']:7.3f} Z={r['Z']:7.3f}  P<= {out[0]:.3f} {out[1]:.3f} {out[2]:.4f}  a(P<0.1)={hi:.3f}  box |Im|<={math.pi/(3*R):.3f}")
# asymptotic: P_bad <~ (log q) * (alpha R/2)^2 /(R e^{alpha R}) with R=2 log q -> q^{-2 alpha} (log q)^2 alpha^2/2
for q in (1e6,1e9,1e12):
    R=2*math.log(q); Zup=math.log(q/math.pi)  # approx (arch ~ log q - const, pole term 2S0/q ~ 4/R)
    print(f"q={q:.0e}: P_bad(0.05)<={Zup/h_i(0.05,R):.3g}  P_bad(0.1)<={Zup/h_i(0.1,R):.3g}  P_bad(0.25)<={Zup/h_i(0.25,R):.3g}")
