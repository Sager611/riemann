# uc_t4.py - explore2 T4: Poisson/Harnack bound for the vertical derivative of S = -Im(H'/H)(x+ih) (positive harmonic above the zeros)
import mpmath as mp, pickle, time
mp.mp.dps = 40
T0 = time.time()
d = pickle.load(open('jets.pkl', 'rb'))
c = [mp.mpf(x) for x in d['cH']][:71]   # t = 0 Taylor coefficients of H_0/S at X0 (z = 2 gamma ~ 14010.163)
I = mp.mpc(0, 1)
def pv(cc, z):
    a = mp.mpc(0)
    for x in reversed(cc): a = a*z + x
    return a
c1 = [(n+1)*c[n+1] for n in range(len(c)-1)]; c2 = [(n+1)*c1[n+1] for n in range(len(c1)-1)]
print("h   max_x h*dS/dh/S   min_x   mean_x   (x in [-3,3], step 0.02; Poisson bound: h*dS/dh <= S since zeros real at t=0)")
for h in [mp.mpf('0.69'), mp.mpf(1), mp.mpf('1.5'), mp.mpf('1.88')]:
    rs = []; Ss = []
    for k in range(301):
        x = mp.mpf(-3) + k*mp.mpf('0.02'); z = x + I*h
        p, p1, p2 = pv(c, z), pv(c1, z), pv(c2, z)
        f = p1/p; fp = (p2*p - p1*p1)/(p*p)
        S = -mp.im(f); dS = -mp.re(fp)
        rs.append(h*dS/S); Ss.append(S)
    print(mp.nstr(h,3), mp.nstr(max(rs),4), mp.nstr(min(rs),4), mp.nstr(sum(rs)/len(rs),4), " S range", mp.nstr(min(Ss),4), mp.nstr(max(Ss),4), " argmax x=", mp.nstr(mp.mpf(-3)+rs.index(max(rs))*mp.mpf('0.02'),3))
print("t=%.1fs" % (time.time()-T0))
