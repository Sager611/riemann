
import mpmath as mp
mp.mp.dps = 50
lg = [mp.mpf(l.split()[1]) for l in open("gammas_log.txt")]
g = [mp.e**x for x in lg]
# sanity: Xi(t) ~ sum (-1)^k g_k t^{2k}
mp.mp.dps=50
for t in [14.0, 14.134725, 14.3, 21.022040]:
    print("Xi-sum", t, mp.nstr(sum((-1)**k*g[k]*mp.mpf(t)**(2*k) for k in range(len(g))),8))
N = len(g)-1
L1 = {n: mp.e**(2*lg[n]-lg[n-1]-lg[n+1]) for n in range(1,N)}
eps = {n: 1-1/L1[n] for n in L1}
C1 = {}
for n in range(2,N-1):
    Th = eps[n]**2/(eps[n-1]*eps[n+1])
    C1[n] = n**2*abs(mp.log(Th))
for n in [10,20,30,40,50,60,80,100,120,140,157]:
    print(n, "a(n)=n*eps", mp.nstr(n*eps[n],8), "Theta1<1?", C1[n]>0, "C1", mp.nstr(C1[n],12), "n*(C1-1)", mp.nstr(n*(C1[n]-1),8))
print("C1(50)-C1(30) =", mp.nstr(C1[50]-C1[30],8), " paper Lemma 11.6 claims <1e-9")
print("C1(150)-C1(100) =", mp.nstr(C1[150]-C1[100],8))
