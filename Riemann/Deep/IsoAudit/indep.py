# indep.py: independent float (mpmath, 30 digits) recomputation of the isolated-zeros numbers. Not a certificate.
# Run: cd Riemann/Deep/IsoAudit && OMP_NUM_THREADS=1 uv run --no-project --with mpmath python3 indep.py
from indep_lib import *
S2 = mp.sqrt(2)
f0, Kr0, Ki0 = kern(S2)
C0 = Cq(f0); CMT = mp.mpf(1)/2 + mp.cot(1/S2)/S2
print('C(f0) quadrature =', mp.nstr(C0, 20), ' closed 1/2+cot(1/sqrt2)/sqrt2 =', mp.nstr(CMT, 20))
for u in ('0.002', '0.3'):
    u = mp.mpf(u); print('K(%s): closed %s quad %s' % (u, mp.nstr(Kr0(u), 18), mp.nstr(2*mp.quad(lambda x: f0(x)*mp.cos(2*PI*u*x), [0, 0.5]), 18)))
for y in ('0.5', '2.3', '5'):
    y = mp.mpf(y); print('K(i%s): closed %s quad %s' % (y, mp.nstr(Ki0(y), 18), mp.nstr(2*mp.quad(lambda x: f0(x)*mp.cosh(2*PI*y*x), [0, 0.5]), 18)))
# --- paper reproduction ---
H0 = Hfun(Ki0, mp.mpf(0))
a0, y0 = supratio(lambda y: H0(y), mp.mpf(0), mp.mpf('0.1'), 5, mp.mpf('0.002'))
print('alpha_0 = max H0(y)/y on [0.1,5] = %s at y = %s (paper < 0.21462; cert [0.214612853119, 0.214612853122])' % (mp.nstr(a0, 15), mp.nstr(y0, 8)))
print('2 - C_MT - 0.0806 alpha_0 =', mp.nstr(2 - CMT - mp.mpf('0.0806')*a0, 12), '(paper > 0.655202)')
print('5/2 - 3/2 C_MT =', mp.nstr(mp.mpf(5)/2 - 3*CMT/2, 12), '(paper 0.50875)')
tau0 = mp.findroot(lambda t: Kr0(2*t) - mp.mpf(5)/6, 0.17)
print('tau_0 (K0(2 tau) = 5/6) =', mp.nstr(tau0, 12), '(paper 0.169641...) > 1/6:', tau0 > mp.mpf(1)/6)
print('H0(2) =', mp.nstr(H0(2), 12), '(paper > 0.4092)')
# --- column A at c = 0.001 (paper inputs, D = 0.0806, our explicit errors) ---
c = mp.mpf('0.001'); lam = lamC(Kr0, c); Hc = Hfun(Ki0, c)
aA, _ = supratio(lambda y: Hc(y) - 3*lam, mp.mpf(0), mp.mpf('0.1'), 8, mp.mpf('0.002'))
print('column A c=0.001: lam =', mp.nstr(lam, 10), ' alpha~ =', mp.nstr(aA, 12), ' p =', mp.nstr((2-lam-CMT-mp.mpf('0.0806')*aA)/(1-lam), 10), '(cert 0.655200)')
# --- hinge at c = 0.001 ---
print('hinge c = 0.001, back = 0.02 (cert: rho 1.2036, alpha 0.670124, B 1.404e-05, p 0.672487):')
p, d = p_hinge(c, C=CMT, verbose=True)
for bk in ('0.001', '0.005', '0.05', '0.1', '0.2'):
    p, d = p_hinge(c, back=mp.mpf(bk), C=CMT); print('  back %s: p = %s' % (bk, mp.nstr(p, 10)))
# --- c* via Conrey's c_1 with P = sinh(bx)/sinh(b) and with a polynomial P ---
R = 2*PI*mp.mpf('1.2036'); b = TH*R
P = lambda x: mp.sinh(b*x)/mp.sinh(b); dP = lambda x: b*mp.cosh(b*x)/mp.sinh(b)
c1 = 1 + (1/TH)*((mp.exp(2*R)-1)/(2*R))*mp.quad(lambda x: (dP(x)+TH*R*P(x))**2, [0, 1])
print('c_+ from c_1 with sinh P:', mp.nstr(mp.exp(-2*R)*c1, 15), ' closed c*:', mp.nstr(cstar(R), 15))
co = [b**k/mp.factorial(k) for k in range(1, 26, 2)]   # odd Taylor polynomial of sinh(bx), renormalised: P(0)=0, P(1)=1
Pp = lambda x: sum(a*x**(2*j+1) for j, a in enumerate(co))/sum(co)
dPp = lambda x: sum(a*(2*j+1)*x**(2*j) for j, a in enumerate(co))/sum(co)
c1p = 1 + (1/TH)*((mp.exp(2*R)-1)/(2*R))*mp.quad(lambda x: (dPp(x)+TH*R*Pp(x))**2, [0, 1])
print('c_+ with degree-25 polynomial P:', mp.nstr(mp.exp(-2*R)*c1p, 15), ' B:', mp.nstr(mp.log(mp.exp(-2*R)*c1p)/(4*PI), 10), 'vs', mp.nstr(Bh(mp.mpf('1.2036')), 10))
print('elapsed %.1f s' % (time.time()-t0)); sys.stdout.flush()
# --- thresholds (back = 0.02, as cert) ---
def cmax(target, lo, hi, steps=22):
    for _ in range(steps):
        m = (lo+hi)/2
        if p_hinge(m, C=CMT)[0] >= target: lo = m
        else: hi = m
    return lo
for target, lo, hi, claim in ((mp.mpf('0.6724'), mp.mpf('0.003'), mp.mpf('0.0065'), '0.00476'), (mp.mpf('0.672'), mp.mpf('0.009'), mp.mpf('0.013'), '0.01106'),
                              (mp.mpf('0.67'), mp.mpf('0.02'), mp.mpf('0.03'), '0.02485'), (mp.mpf('0.66'), mp.mpf('0.05'), mp.mpf('0.06'), '0.05495'), (mp.mpf('0.6552'), mp.mpf('0.06'), mp.mpf('0.07'), '0.06427')):
    cm = cmax(target, lo, hi)
    pc = p_hinge(mp.mpf(claim), C=CMT)[0]
    print('target %s: c_max = %s (claim %s); p(claim) = %s >= target: %s' % (target, mp.nstr(cm, 8), claim, mp.nstr(pc, 10), pc >= target)); sys.stdout.flush()
# --- c = 1/6 ---
c6 = mp.mpf(1)/6
print('c = 1/6, f_0:'); p_hinge(c6, C=CMT, verbose=True)
w = mp.sqrt(mp.mpf('3.75')); fw, Krw, Kiw = kern(w)
print('c = 1/6, w^2 = 3.75: 2K(c)^2 =', mp.nstr(2*Krw(c6)**2, 8), ' 4c =', mp.nstr(4*c6, 6), ' C(f_w) =', mp.nstr(Cq(fw), 12))
p_hinge(c6, w=w, verbose=True)
print('elapsed %.1f s' % (time.time()-t0))
