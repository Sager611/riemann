"""make_config.py -- parameters of an oocert run (weilooc, 2026-11-05).
Given a tag and a weight file: s = dyadic just above the rigorous Schur bound (bc.schur), gamma = s + gap rounded up to
1/256, Omega = smallest integer with C = a(Omega) - gamma >= Ctarget, u_j = floor(C/E_j) - 1, N = even ~ ratio * Omega L.
usage: make_config.py TAG WFILE OUTDIR [ratio=1.31] [tau_exp=400] [gap=0.8579] [Ctarget=0.1154]"""
import sys, os, math, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from flint import arb, ctx, fmpq
import bc
tag, wfile, out = sys.argv[1], sys.argv[2], sys.argv[3]
ratio = float(sys.argv[4]) if len(sys.argv) > 4 else 1.31
tau_exp = int(sys.argv[5]) if len(sys.argv) > 5 else 400
gap = float(sys.argv[6]) if len(sys.argv) > 6 else 0.8579
Ct = float(sys.argv[7]) if len(sys.argv) > 7 else 0.1154
ctx.prec = 128
S = bc.Setup(tag); L = S.L_()
bound = bc.schur(S, wfile)[0]
s_num = math.ceil(float(bound.upper())*2**20) + 1
s = fmpq(s_num, 2**20)
g_num = math.ceil((float(s_num)/2**20 + gap)*256)
gamma = arb(fmpq(g_num, 256))
lo, hi = 100, 100000
while hi - lo > 1:
    mid = (lo + hi)//2
    if bc.a_fun(arb(mid)) - gamma >= Ct: hi = mid
    else: lo = mid
Om = hi; C = bc.a_fun(arb(Om)) - gamma; c = arb(Om)*L
E0 = (c + 1)/(arb.pi()*c**2); E1 = (3*c**2 + 6*c + 2)/(arb.pi()*c**3)
u0 = int(float((C/E0).lower())) - 1; u1 = int(float((C/E1).lower())) - 1
N = 2*round(ratio*float(c.mid())/2)
cfg = {'tag': tag, 'N': N, 'F': 448, 'T': 480, 'FV': 512, 'Omega': str(Om), 'gamma': f"{g_num}/256", 's': f"{s_num}/1048576",
       'tau': f"2^-{tau_exp}", 'u0': u0, 'u1': u1, 'wfile': os.path.abspath(wfile), 'block': 480, 'budget': 2100,
       'band_target_bits': 440, 'X_tolbits': 240}
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'config.json'), 'w') as fh: json.dump(cfg, fh, indent=1)
print(f"{tag}: L={float(L.mid()):.6f} Schur {bound.str(12)} -> s={s_num}/2^20, gamma={g_num}/256 (gap {g_num/256 - s_num/2**20:.5f}), Omega={Om}, C={float(C.mid()):.5f}, c={float(c.mid()):.1f}, u=({u0},{u1}), N={N} (N/c {N/float(c.mid()):.3f})")
