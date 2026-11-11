# Adversarial small sets E for Prop SmallSetZeros on real zeros: remove the delta*N zeros with the largest local
# kernel mass c_i = sum_j |g0^((g_i-g_j)L)| w_ij, compare with random E; Prop claims |S(Z-E)-S(Z)|/N -> 0 as delta -> 0.
import sys, gzip, math
import numpy as np
from num_pairs import g0hat, pair_sum, TWO_PI, load
gam_all = load(sys.argv[1])
N = int(sys.argv[2]); gam = gam_all[:N]; T = gam[-1] + 1e-9; L = math.log(T)/TWO_PI
c = np.full(N, g0hat(np.array([0.0]))[0]); d = 1
while True:
    dl = gam[d:] - gam[:-d]
    if dl.min() > 100: break
    k = np.abs(g0hat(dl*L))*4/(4 + dl*dl); c[d:] += k; c[:-d] += k; d += 1
S, _ = pair_sum(gam, L, g0hat, 100.0)
rng = np.random.default_rng(0)
print(f'T={T:.1f} N={N}: S(Z)/N={S/N:.5f}; local mass c_i: mean={c.mean():.3f} max={c.max():.3f}', flush=True)
for delta in (1e-3, 1e-2, 1e-1):
    m = int(delta*N)
    top = np.argsort(c)[-m:]; keep = np.ones(N, bool); keep[top] = False
    Sa, _ = pair_sum(gam[keep], L, g0hat, 100.0)
    keep = np.ones(N, bool); keep[rng.choice(N, m, replace=False)] = False
    Sr, _ = pair_sum(gam[keep], L, g0hat, 100.0)
    print(f'  delta={delta:g}: adversarial (S(Z-E)-S(Z))/N={(Sa-S)/N:+.5f} (/delta={(Sa-S)/N/delta:+.3f});  random: {(Sr-S)/N:+.5f} (/delta={(Sr-S)/N/delta:+.3f})', flush=True)
