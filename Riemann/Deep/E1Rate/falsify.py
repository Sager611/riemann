# e1rate falsification: (F1) Mellin separation cost of the Bessel phase; (F2) Lemma I (incomplete Kloosterman sum with Bessel phase).
import numpy as np, math, sys
rng = np.random.default_rng(1)
def bump(z):  # smooth bump on (1,2)
    out = np.zeros_like(z); m = (z > 1) & (z < 2); s = (z[m]-1)*(2-z[m]); out[m] = np.exp(-1/(4*s)) ; return out
print("F1: Phi(x)=exp(i*X0*sqrt(x))*bump(x) on x in (1,2); L1 norm of Mellin transform over t vs sup|Phi|; predicted ratio ~ sqrt(X0)")
for X0 in [0, 10, 40, 160, 640, 2560]:
    M = 2**17; u = np.linspace(0, math.log(2.0), M, endpoint=False); du = u[1]-u[0]
    x = np.exp(u); f = np.exp(1j*X0*np.sqrt(x))*bump(x)
    # Mellin at s = it: int f(e^u) e^{itu} du ; FFT over a long zero-padded grid
    L = 2**20; F = np.fft.fft(np.concatenate([f, np.zeros(L-M)]))*du
    dt = 2*math.pi/(L*du); l1 = np.sum(np.abs(F))*dt/(2*math.pi)
    print(f"  X0={X0:6d}  int|Mellin|dt/2pi = {l1:.5f}  sup|Phi| = {np.max(np.abs(f)):.5f}  ratio/sqrt(1+X0) = {l1/np.max(np.abs(f))/math.sqrt(1+X0):.4f}", flush=True)
def isprime(n): return n > 1 and all(n % d for d in range(2, int(n**.5)+1))
print("F2: S = sum_k e(c*kbar/p) w(k/D) exp(i*X0*sqrt(D/k)); Lemma I bound B = D/p + sqrt(p*(1+X0)); trivial = D")
for p in [1009, 4001]:
    inv = None
    for D in [int(p**0.75), p, 4*p]:
        k = np.arange(D, 2*D+1); kk = k[k % p != 0]
        kbar = np.array([pow(int(t), -1, p) for t in kk]); w = bump(kk/D)
        for X0 in [0, 30, 300, 3000]:
            vals = []
            for c in rng.integers(1, p, size=12):
                ph = 2*math.pi*((int(c)*kbar) % p)/p + X0*np.sqrt(D/kk)
                vals.append(abs(np.sum(w*np.exp(1j*ph))))
            B = D/p + math.sqrt(p*(1+X0))
            print(f"  p={p} D={D} X0={X0:5d}  mean|S|={np.mean(vals):9.2f} max|S|={np.max(vals):9.2f}  B={B:9.2f}  max/B={np.max(vals)/B:.3f}  trivial={np.sum(w):9.1f}", flush=True)
