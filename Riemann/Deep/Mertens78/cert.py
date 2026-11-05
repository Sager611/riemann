# cert.py -- rigorous (python-flint arb) certificate checker for the mertens78 bound.
# Usage: python3 cert.py cert.json
# cert.json: {"sigma1": "0.92", "eta": "0.02", "masses": [["1.01","0","3.0"], ...] (sigma_j, y_j, lambda_j exact decimals),
#             "T1": "1e4", "harnack": true/false, "xa": "1e100", "ch": "1.0", "theta": "0.93"}
import json, sys
from flint import arb, acb, ctx
ctx.prec = 120
B = arb(7)/8
PI = arb.pi()

def A_(s): return arb(s)

class Cert:
    def __init__(self, d):
        self.d = d
        self.s1 = A_(d["sigma1"]); self.eta = A_(d["eta"])
        if d.get("harnack"):
            s0 = A_(d["sigma0"])
            self.m = [(s0, arb(0), (s0 - B)/(self.s1 - B))]   # Harnack: lambda = (x0/x1) delta_0
        else:
            self.m = [(A_(a), A_(b), A_(c)) for a, b, c in d["masses"]]
        self.kappa = (arb(1)/8 + self.eta)/2
        self.A = (3*(1+self.eta).zeta()).log() + self.kappa*(arb(8)/15).log()
        self.Lam = sum((lam for _, _, lam in self.m), arb(0))
        self.a = (self.Lam - 1)*self.kappa
        self.x1 = self.s1 - B
    # harmonic majorant h(s), s = sig + i t
    def h(self, sig, t):
        return self.A + (1+self.kappa)*((1+sig)**2 + t**2).log()/2 - (sig**2 + t**2).log()/2
    def L0(self, sig, t):
        return (sig.zeta()/(2*sig).zeta()).log() + (sig**2 + t**2).log()/2 - ((sig-1)**2 + t**2).log()/2
    def L1(self, sig, t):
        return sig.zeta().log() + ((sig-1)**2 + t**2).log()/2 - (sig**2 + t**2).log()/2
    def logPhi(self, t):
        s1 = self.s1
        tot = -self.h(s1, t) + ((s1-1)**2 + t**2).log()/2 - (s1**2 + t**2).log()/2
        for sig, y, lam in self.m:
            if lam > 0:
                tot += lam*(self.h(sig, t+y) + self.L0(sig, t+y))
            else:
                tot += -lam*(-self.h(sig, t+y) + self.L1(sig, t+y))
        return tot
    # ---- condition (C1): sum lam_j P_{x_j}(eta - y_j) >= P_{x1}(eta) for all real eta; (C2): sum lam_j x_j >= x1
    def F(self, e):
        tot = -self.x1/(e**2 + self.x1**2)
        for sig, y, lam in self.m:
            x = sig - B
            tot += lam*x/((e-y)**2 + x**2)
        return tot
    def check_C(self, Y=None, maxdepth=40):
        x1 = self.x1
        c2 = sum((lam*(sig-B) for sig, y, lam in self.m), arb(0)) - x1
        assert c2 > 0, "C2 fails"
        cmax = max(abs(float(y.mid())) for _, y, _ in self.m)
        if self.d.get("harnack"):
            # single mass rho*delta_0 at x0 with rho = x0/x1: analytic (see PROOF.md)
            return "harnack (analytic)", 0
        def bracket(Yv):
            Ya = arb(Yv); br = -x1
            for sig, y, lam in self.m:
                x = sig - B; c = abs(y)
                if lam > 0: br += lam*x*Ya**2/((Ya + c)**2 + x**2)
                else: br += lam*x*Ya**2/((Ya - c)**2)
            return br
        # tail bracket at |eta| = Y (nondecreasing in |eta| > max c_j; see PROOF.md)
        Y = None
        for Yc in [2*cmax + 50, 1e3, 1e4, 1e5, 1e6, 1e7]:
            if Yc > cmax and bracket(Yc) > 0: Y = Yc; break
        assert Y is not None, "tail bracket fails"
        self.Ytail = Y
        # adaptive bisection on [-Y, Y]: uniform on [-Y0, Y0], geometric (ratio 1.002) beyond
        Y0 = min(Y, 2*cmax + 50)
        stack = [(-Y0 + 2*Y0*k/4000, -Y0 + 2*Y0*(k+1)/4000, 0) for k in range(4000)]
        w = Y0
        while w < Y:
            w2 = min(w*1.002, Y); stack.append((w, w2, 0)); stack.append((-w2, -w, 0)); w = w2
        n = 0; worst = None
        while stack:
            lo, hi, dep = stack.pop(); n += 1
            iv = arb((lo+hi)/2, (hi-lo)/2)
            v = self.F(iv)
            if v > 0:
                continue
            if dep >= maxdepth:
                raise AssertionError("C1 not certified near %r" % lo)
            mid = (lo+hi)/2
            stack.append((lo, mid, dep+1)); stack.append((mid, hi, dep+1))
        return "bisection ok", n
    # ---- tail constants for t >= T1
    def tail(self, T1):
        T1 = arb(T1)
        R = self.s1
        K = (self.Lam - 1)*self.A
        cmax = 0.0
        for sig, y, lam in self.m:
            c = abs(y); cmax = max(cmax, float(c.upper()))
            if lam > 0:
                R += lam*((1+self.kappa)*(c + 1 + sig) + 2*c + 2)
                K += lam*(sig.zeta()/(2*sig).zeta()).log()
            else:
                R += (-lam)*((3 + 2*self.kappa)*c + sig)
                K += (-lam)*sig.zeta().log()
        assert T1 >= 2*cmax and T1 >= 1
        return K + R/T1   # log C_inf
    # ---- J0 = 2 int_0^T1 Phi(t)/|s1| dt (upper bound)
    def J0(self, T1, w0=0.002, rel=0.002):
        T1f = float(T1); t = 0.0; tot = arb(0); n = 0
        while t < T1f:
            w = max(w0, rel*t); t2 = min(t + w, T1f)
            iv = arb((t+t2)/2, (t2-t)/2)
            val = (self.logPhi(iv).exp()/((self.s1**2 + iv**2).sqrt()))
            tot += arb(val.upper())*(arb(t2) - arb(t)); t = t2; n += 1
        return 2*tot, n

def J0_direct(s1, T1):
    """Rigorous upper bound for 2*int_0^T1 |1/zeta(s1+it)|/|s1+it| dt via arb ball evaluations of zeta."""
    s1 = arb(s1); T1f = float(T1); t = 0.0; tot = arb(0); n = 0; mx = 0.0
    while t < T1f:
        w = 0.02 if t < 100 else 0.01
        t2 = min(t + w, T1f)
        iv = arb((t+t2)/2, (t2-t)/2)
        zv = acb(s1, iv).zeta()
        lo = zv.abs_lower()
        assert lo > 0
        inv = 1/lo
        mx = max(mx, float(inv.upper()))
        tot += arb((inv/((s1**2 + arb(t)**2).sqrt())).upper())*(arb(t2) - arb(t)); t = t2; n += 1
    return 2*tot, n, mx

def final(c, T1, xa, ch, theta, J0=None):
    """Check |M(x)| <= C x^theta for all x >= xa; returns C (upper bound)."""
    a = c.a; s1 = c.s1
    theta = arb(theta); xa = arb(xa); ch = arb(ch); T1a = arb(T1)
    logCinf = c.tail(T1)
    Cinf = logCinf.exp()
    if J0 is None:
        J0, n = c.J0(T1)
    else:
        n = -1
    # exponent check: sigma1 + a(1-theta) - theta <= 0
    expo = s1 + a*(1-theta) - theta
    assert expo <= 0, "theta below theta*"
    h = ch*xa**theta
    assert h <= xa
    X = xa + h/2; tau = 2*X/h
    # J(tau) <= J0 + 2 C_inf (1/a + 1/(1-a)) tau^a for every tau > 0 (PROOF.md, Lemma B5)
    T1t = ch/4
    T2t = xa**(-theta)
    J0eff = J0
    T3t = (X**s1/(2*PI*xa**theta))*(J0eff if J0eff > 0 else arb(0))
    T4t = (1 + (ch/2)*xa**(theta-1))**s1 * ((2/ch) + xa**(theta-1))**a * xa**expo * 2*Cinf*(1/a + 1/(1-a))/(2*PI)
    C = T1t + T2t + T3t + T4t
    return dict(C=C, J0=J0, nJ0=n, logCinf=logCinf, a=a, T1t=T1t, T2t=T2t, T3t=T3t, T4t=T4t, J0eff=J0eff)

if __name__ == "__main__":
    d = json.load(open(sys.argv[1]))
    c = Cert(d)
    print("kappa", c.kappa, "A", c.A, "Lam", c.Lam, "a", c.a)
    print("C1/C2:", c.check_C())
    r = final(c, d["T1"], d["xa"], d["ch"], d["theta"])
    for k, v in r.items(): print(k, v)

