# lambda_n = n [z^n] log xi(1/(1-z)) by DFT on |z| = r (ball arithmetic, python-flint).
# usage: python3 fft_lambda.py r M prec nmax out.txt [gamma]
# mode gamma: Maslanka's trend, n [z^n] log(pi^{-s/2} Gamma(1+s/2)) (no zeta; entire in Re s > 0)
# Unconditional for the computed range: the image of |z| <= r is a disk in Re s > 1/(1+r) > 1/2 with
# |Im s| <= r/(1-r^2); xi has no zeros there (RH verified far above that height).
# Aliasing error |sum_{j>=1} a_{n+jM} r^{jM}| <= r^M * sum ... is printed as a bound (a_k <= log k for k >= 2).
import sys, time, math
from flint import acb, arb, ctx

r_s, M, prec, nmax, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
GAM = len(sys.argv) > 6 and sys.argv[6] == 'gamma'
ctx.prec = prec
r = arb(r_s)
pi = arb.pi(); logpi = pi.log(); log2 = arb(2).log(); twopi = 2 * pi
t0 = time.time()
half = M // 2
vals = [None] * M
# Branches: s, s-1 stay in the upper half plane on the upper arc (principal logs continuous); lgamma is the
# analytic log Gamma; log zeta: principal where Re s >= 2 (Re zeta > 0.355 there), else continued from the
# previous point (on that part of the arc |ds| per step is small and no zero of zeta is within distance ~1).
prev = None; maxjump = 0.0; nunwrap = 0
def xi_log_real(s):  # s real > 1/2: xi(s) > 0
    if GAM: return -s * logpi / 2 + (1 + s / 2).lgamma()
    return (s.log() + abs(s - 1).log() - log2 - s * logpi / 2 + (s / 2).lgamma() + abs(s.zeta()).log())
for k in range(1, half):
    ang = twopi * k / M
    z = acb(r * ang.cos(), r * ang.sin())
    s = 1 / (1 - z)
    if GAM:
        vals[k] = -s * logpi / 2 + (1 + s / 2).lgamma(); continue
    lz = s.zeta().log()
    if float(s.real.mid()) < 2:
        d = float((lz.imag - prev.imag).mid())
        kk = round(d / (2 * math.pi))
        if kk:
            lz = lz - acb(0, twopi * kk)
        maxjump = max(maxjump, abs(float((lz.imag - prev.imag).mid())))
        nunwrap += 1
    prev = lz
    vals[k] = s.log() + (s - 1).log() - log2 - s * logpi / 2 + (s / 2).lgamma() + lz
    if k % 16384 == 0:
        print(f'k={k} t={time.time()-t0:.1f}s', flush=True)
vals[0] = acb(xi_log_real(1 / (1 - r)))
vals[half] = acb(xi_log_real(1 / (1 + r)))
# continuity check at the endpoint: last interior point vs exact real value
end_im = float(vals[half - 1].imag.mid())
print(f'Im f at last interior point {end_im:.3e}; Re jump {float((vals[half].real - vals[half-1].real).mid()):.3e}', flush=True)
end_im = end_im if abs(end_im) > 1e-2 else 0.0
print(f'unwrapped points {nunwrap}; max |d arg zeta| there {maxjump:.3e}; Im f(-r) = {end_im:.3e} (must be 0)', flush=True)
assert abs(end_im) < 1e-6 and maxjump < 0.5
vals[half] = acb(vals[half].real, 0)
vals[0] = acb(vals[0].real, 0)
for k in range(1, half):
    vals[M - k] = vals[k].conjugate()
c = acb.dft(vals)
print(f'dft done t={time.time()-t0:.1f}s', flush=True)
alias = float((r ** M).mid()) * 4 * math.log(M + nmax)  # aliasing in a_n: |sum_j a_{n+jM} r^{jM}| (a_k <= log k)
with open(out, 'w') as fh:
    fh.write(f'# r={r_s} M={M} prec={prec}\n# n lambda_n radius alias_bound\n')
    rn = arb(1)
    rinv = 1 / r
    worst = 0.0
    for n in range(0, nmax + 1):
        a = c[n].real / M * rn
        lam = a * n
        rad = float(lam.rad())
        ab = alias * n
        worst = max(worst, rad + ab)
        fh.write(f'{n} {lam.mid().str(30, radius=False)} {rad:.2e} {ab:.2e}\n')
        rn = rn * rinv
print(f'wrote {out}; worst error bound {worst:.3e}; total {time.time()-t0:.1f}s', flush=True)
