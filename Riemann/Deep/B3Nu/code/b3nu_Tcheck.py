import sys, math, time, numpy as np
from b3nu_T import *
t0 = time.time()
for T, th in [(1e6, 0.45), (1e6, 0.40), (1e6, 0.47)]:
    H = T ** th; tau = T / H ** 2
    TI = Tinterp(T, H, 1.0 / 2000, 4.5 * tau)
    rng = np.random.default_rng(1); xs = np.exp(rng.uniform(math.log(1/1500), math.log(4 * tau), 12))
    w, Phi = TI.full; ws, Ps = TI.sub
    ex = Tbr(T, H, xs, w, Phi, 1) + Tbr(T, H, xs, w, Phi, -1)
    sb = Tbr(T, H, xs, ws, Ps, 1) + Tbr(T, H, xs, ws, Ps, -1)
    it = TI(xs); scale = H * (8 * np.pi * xs * T) ** -0.25   # H X0^{-1/2}, X0 = 2 sqrt(2 pi) sqrt(T x)
    print(f"T={T:.0e} th={th} H={H:.1f} tau={tau:.2f} nw={len(w)} nsub={len(ws)}: max|sub-full|/scale={np.max(np.abs(sb-ex)/scale):.1e} "
          f"max|interp-full|/scale={np.max(np.abs(it-ex)/scale):.1e}  max|T|/scale={np.max(np.abs(ex)/scale):.3f}  t={time.time()-t0:.0f}s", flush=True)
