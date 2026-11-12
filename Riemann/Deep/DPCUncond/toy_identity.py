#!/usr/bin/env python3
"""Toy test of the hairpin identity for DPC's weighted moments with planted off-line zeros.

F(s) = prod_{off-line rho}(s-rho)(s-1+conj(rho)) * prod_k (s-1/2-i g_k) * zeta(s+1)
log F by horizontal continuation from Re s=+inf: principal Log of each linear factor, principal log zeta(s+1) (|log zeta|<1 on Re s>=3/2).
V(t) = (log F(1/2+i(t+h/2)) - log F(1/2+i(t-h/2)))/pi ; R(t)=sum c_q q^{-it}, Rs(t)=sum conj(c_q) q^{it}; g Gaussian.
Identity (Cauchy on [a,b]x[-1,0] minus vertical slits at x0=gamma-+h/2, depth d=beta-1/2):
 int_a^b F dx = int_{a-i}^{b-i} F dt + int_{b-i}^{b} F dt - int_{a-i}^{a} F dt - i sum_slits int_0^d (F_R-F_L)(x0-iy) dy,
 F = g V^j R Rs,  F_R - F_L = g R Rs [(V_L + J)^j - V_L^j],  J = +2i at gamma-h/2, -2i at gamma+h/2.
"""
import sys, json, time
import mpmath as mp
mp.mp.dps = 30

def mob(n):
    m, p, r = n, 2, 1
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    if m > 1: r = -r
    return r

def setup(Q, offl, onl, h, tc, Hw):
    cs = [(q, mob(q)/mp.sqrt(q)) for q in range(1, Q+1) if mob(q) != 0]
    zs = []
    for (d, gam) in offl:
        zs.append(mp.mpc(0.5+d, gam)); zs.append(mp.mpc(0.5-d, gam))
    for gam in onl:
        zs.append(mp.mpc(0.5, gam))
    def logF(s):
        v = mp.log(mp.zeta(s+1))
        for z in zs:
            v += mp.log(s - z)
        return v
    def V(t):
        return (logF(mp.mpf(0.5) + 1j*(t + h/2)) - logF(mp.mpf(0.5) + 1j*(t - h/2)))/mp.pi
    def R(t):
        return mp.fsum(c*mp.power(q, -1j*t) for q, c in cs)
    def Rs(t):
        return mp.fsum(c*mp.power(q, 1j*t) for q, c in cs)
    def g(t):
        return mp.exp(-(t-tc)**2/(2*Hw**2))
    return cs, V, R, Rs, g

def run(Q, offl, onl, h, tc, Hw, js, a, b):
    cs, V, R, Rs, g = setup(Q, offl, onl, h, tc, Hw)
    brk = sorted(set([a, b] + [gg + s*h/2 for gg in onl for s in (-1, 1)] + [gg + s*h/2 for (_, gg) in offl for s in (-1, 1)]))
    brk = [x for x in brk if a <= x <= b]
    out = {}
    for j in js:
        Fr = lambda x: g(x)*V(mp.mpf(x))**j*R(x)*Rs(x)
        I_real = mp.quad(Fr, brk, maxdegree=8)
        Fl = lambda x: g(x - 1j)*V(x - 1j)**j*R(x - 1j)*Rs(x - 1j)
        I_low = mp.quad(Fl, [a, b])
        sR = mp.quad(lambda y: g(b - 1j*y)*V(b - 1j*y)**j*R(b - 1j*y)*Rs(b - 1j*y)*(-1j), [1, 0])
        sL = mp.quad(lambda y: g(a - 1j*y)*V(a - 1j*y)**j*R(a - 1j*y)*Rs(a - 1j*y)*(-1j), [1, 0])
        hp = mp.mpc(0)
        eps = mp.mpf(10)**(-22)
        hps = []
        for (d, gam) in offl:
            for (x0, J) in ((gam - h/2, 2j), (gam + h/2, -2j)):
                def jump(y, x0=x0, J=J):
                    t = x0 - 1j*y
                    VL = V(x0 - eps - 1j*y)
                    VR = V(x0 + eps - 1j*y)
                    return g(t)*R(t)*Rs(t)*(VR**j - VL**j), VR - VL
                f = lambda y: jump(y)[0]
                val = -1j*mp.quad(f, [0, d])
                jchk = jump(d/3)[1]
                hps.append((float(x0), complex(val), complex(jchk)))
                hp += val
        rhs = I_low + sR - sL + hp
        out[j] = dict(I_real=complex(I_real), I_low=complex(I_low), sides=complex(sR - sL), hairpins=complex(hp),
                      rhs=complex(rhs), resid=float(abs(I_real - rhs)), resid_nohp=float(abs(I_real - (I_low + sR - sL))),
                      slit_detail=hps)
    return out

if __name__ == "__main__":
    t0 = time.time()
    h = mp.mpf("0.8"); tc = mp.mpf(50); Hw = mp.mpf(3)
    a, b = tc - 30, tc + 30
    onl = [mp.mpf(x) for x in ("47.0", "49.5", "52.0", "53.2")]
    cases = {
      "single": [(mp.mpf("0.25"), mp.mpf("51.3"))],
      "cluster": [(mp.mpf("0.25"), mp.mpf("51.3")), (mp.mpf("0.1"), mp.mpf("51.65")), (mp.mpf("0.4"), mp.mpf("48.2"))],
    }
    res = {}
    for name, offl in cases.items():
        for Q in (1, 12):
            r = run(Q, offl, onl, h, tc, Hw, (1, 2, 3), a, b)
            res[name + "_Q" + str(Q)] = r
            for j, v in r.items():
                print(name, "Q=%d j=%d" % (Q, j), "I_real=%.10g%+.10gi" % (v["I_real"].real, v["I_real"].imag),
                      "hairpins=%.6g%+.6gi" % (v["hairpins"].real, v["hairpins"].imag),
                      "resid=%.2e" % v["resid"], "resid_without_hairpins=%.2e" % v["resid_nohp"], flush=True)
    def enc(o):
        if isinstance(o, complex): return [o.real, o.imag]
        if isinstance(o, dict): return {str(k): enc(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [enc(v) for v in o]
        return o
    json.dump(enc(res), open(sys.argv[1] if len(sys.argv) > 1 else "toy_identity.json", "w"), indent=1)
    print("time %.1f s" % (time.time() - t0))
