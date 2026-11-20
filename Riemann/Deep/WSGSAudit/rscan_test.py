# rscan_test.py (wsgsaudit): empirical test of rscan.c's a-priori rounding bound E (cert_simple.py l. 115-170):
# runs rscan_dbg.c (rscan.c + printing the prime sum s_c(k)) and compares 2|s_c(k) - s(k)| with E, s(k) = sum_n w_n cos(r_k ln n) in Arb.
import sys, math, subprocess
sys.path.insert(0, '/tmp/wsgs_copy')
import mpmath as mpm
from flint import arb, ctx, fmpq
import cert_simple as CS
cq = fmpq(int(sys.argv[1])); ctx.prec = 192
S2 = CS.SymZeta(cq); pi = arb.pi(); L = arb(cq).log(); two64 = arb(2)**64; u = 2.0**-53
tab = [f"{float((2*pi*a/65536).cos().mid()).hex()} {float((2*pi*a/65536).sin().mid()).hex()}" for a in range(65536)]
K = float((2*pi/two64).mid()); mpm.mp.dps = 60; Lm = mpm.log(mpm.mpf(int(cq.p))/mpm.mpf(int(cq.q)))
at = []; emax = arb(0); dW = arb(0); sigW = 0.0
for n, _, w in S2.atoms:
    Bn = int(mpm.floor(mpm.log(n)/(2*Lm)*mpm.mpf(2)**64 + mpm.mpf(1)/2)); en = abs(arb(n).log()/(2*L)*two64 - Bn)
    if en.upper() > emax.upper(): emax = arb(en.upper())
    Wn = float(w.mid()); dW += abs(arb(Wn) - w).upper(); sigW += abs(Wn); at.append(f"{Bn % 2**64} {Wn.hex()}")
blocks = [(1, 4001), (10**7 + 1, 10**7 + 2001), (2*10**9 + 1, 2*10**9 + 401), (2, 2002)]
kmax = max(b[1] for b in blocks); na = len(S2.atoms); gam = na*u/(1 - na*u)
E = 2*(sigW*(1 + u)*(gam*(1 + u) + u + 1e-15) + float(dW.upper())*1.0001 + sigW*1.0001*2*math.pi*kmax*float(emax.upper())*2.0**-64)
E += 3*u*(30 + 2*sigW*1.0001); E *= 1.01   # |H|, |T| <= 30 here (H, T are set to 0 below; the term only enlarges E)
inp = "\n".join(tab) + f"\n{K.hex()}\n{na}\n" + "\n".join(at) + "\n" + f"{(1.0).hex()} {(1e-6).hex()}\n{len(blocks)}\n" + "\n".join(f"{a} {b} {(0.0).hex()} {(0.0).hex()}" for a, b in blocks) + "\n"
out = subprocess.run(['/tmp/wsgs_rscan_dbg'], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
worst = 0.0; nn = 0
for line in out:
    if not line: continue
    k, s = line.split(); k = int(k); r = S2.r(k)
    sa = sum((w*(r*ln).cos() for _, ln, w in S2.atoms), arb(0))
    d = 2*abs(arb(float(s)) - sa).upper(); worst = max(worst, float(d)); nn += 1
print(f"c={sys.argv[1]} atoms={na} nodes={nn} kmax={kmax} max 2|s_c - s| = {worst:.3e}  a-priori E = {E:.3e}  ratio {worst/E:.3f}")
