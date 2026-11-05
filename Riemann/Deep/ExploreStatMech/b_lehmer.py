
# B: Lee-Yang lift test of the Riemann-Siegel / AFE main term at Lehmer's pair t~7005.06
import mpmath as mp, sys
mp.mp.dps = 30
def ZN_main(t, N=None):
    t = mp.mpf(t)
    if N is None: N = int(mp.floor(mp.sqrt(t/(2*mp.pi))))
    th = mp.siegeltheta(t)
    return 2*mp.fsum(mp.cos(th - t*mp.log(n))/mp.sqrt(n) for n in range(1, N+1)), N
def F(s, N):
    chi = mp.power(2, s)*mp.power(mp.pi, s-1)*mp.sin(mp.pi*s/2)*mp.gamma(1-s)
    return mp.fsum(mp.power(n, -s) for n in range(1, N+1)) + chi*mp.fsum(mp.power(n, s-1) for n in range(1, N+1))
a, b, h = mp.mpf('7004.9'), mp.mpf('7005.3'), mp.mpf('0.002')
ts = [a + k*h for k in range(int((b-a)/h)+1)]
zt = [mp.siegelz(t) for t in ts]
zm = [ZN_main(t)[0] for t in ts]
def signch(v, ts):
    return [float((ts[i]+ts[i+1])/2) for i in range(len(v)-1) if v[i]*v[i+1] < 0]
print('N =', ZN_main(ts[0])[1])
print('true Z sign changes:', signch(zt, ts))
print('main-term sign changes:', signch(zm, ts))
i0 = min(range(len(ts)), key=lambda i: abs(ts[i]-mp.mpf('7005.08')))
print('Z(7005.08)=', mp.nstr(zt[i0], 6), ' Zmain(7005.08)=', mp.nstr(zm[i0], 6))
# local extremum of Zmain in the pair window
w = [i for i in range(len(ts)) if mp.mpf('7005.0') <= ts[i] <= mp.mpf('7005.15')]
imx = max(w, key=lambda i: zm[i]) ; imn = min(w, key=lambda i: zm[i])
print('Zmain max/min in [7005.0,7005.15]:', mp.nstr(zm[imx],6), 'at', mp.nstr(ts[imx],8), '|', mp.nstr(zm[imn],6), 'at', mp.nstr(ts[imn],8))
imx2 = max(w, key=lambda i: zt[i]); imn2 = min(w, key=lambda i: zt[i])
print('Z max/min in [7005.0,7005.15]:', mp.nstr(zt[imx2],6), 'at', mp.nstr(ts[imx2],8), '|', mp.nstr(zt[imn2],6), 'at', mp.nstr(ts[imn2],8))
N = 33
found = []
for d0 in ['0.001','0.003','0.01','0.03','0.06','0.1']:
    for t0 in ['7005.04','7005.06','7005.08','7005.10','7005.12']:
        try:
            r = mp.findroot(lambda s: F(s, N), mp.mpc(0.5+float(d0), float(t0)))
            if abs(F(r,N)) < mp.mpf('1e-15') and abs(r.imag-7005.08) < 0.5:
                if all(abs(r-q) > 1e-8 for q in found): found.append(r)
        except Exception as e:
            pass
for r in found:
    print('zero of F_33: sigma-1/2 = %s, t = %s, |F| = %s' % (mp.nstr(r.real-0.5, 10), mp.nstr(r.imag, 12), mp.nstr(abs(F(r,N)), 3)))
