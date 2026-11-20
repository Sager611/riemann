# b3audit F4: Lemma R off w = 1. M(s) = int chi(x/T) A_{-1}(x) x^{s-1} dx vs M_lead = (H/(2 sqrt pi)) (T/2pi)^{w-1}, w = s - beta,
# with A_{-1} from the direct y-trapezoid of b3a_lemB.py (chi~ = 1 on the computed range; the profile is < 1e-13 of its max at the ends).
# Claim under test: |M/M_lead - 1| <= C_R H/T uniformly on |w - 1| <= 1/4.
import sys, math, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from b3a_lemB import profile
for T, th in [(1e4, 0.40), (1e5, 0.40), (1e6, 0.40), (1e6, 0.45), (1e6, 0.33)]:
    H = T**th; K = T/H; yb = 6.8; dlo, dhi = math.exp(-2*yb/H) - 1, math.exp(2*yb/H) - 1
    hd = 2*math.pi/(T*max(abs(dlo), dhi)/(1 - abs(dlo)) + 12*K); nd = int((dhi - dlo)/hd) + 2
    dl = np.linspace(dlo, dhi, nd); xi = T*(1+dl)/(2*math.pi); dxi = T*(dl[1]-dl[0])/(2*math.pi)
    for be in [0.0, -0.03+0.07j]:
        A, _ = profile(T, H, be, dl); out = []
        for w in [1.2+0.1j, 0.8-0.15j, 1.0+0.25j, 1.25]:
            M = np.sum(A*xi**(w + be - 1))*dxi; Ml = H/(2*math.sqrt(math.pi))*(T/(2*math.pi))**(w - 1)
            out.append(abs(M/Ml - 1))
        print(f"T={T:.0e} th={th} beta={be}: max_w |M/M_lead-1| = {max(out):.3e}  (H/T = {H/T:.2e}, ratio to H/T = {max(out)/(H/T):.3f}); per w {['%.2e' % o for o in out]}", flush=True)
