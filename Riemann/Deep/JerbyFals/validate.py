# Validate gram_win (double-double phases) and the old plain-double method against mpmath siegelz at 20 random Gram points per decade.
import random, subprocess, sys, time, mpmath
mpmath.mp.dps = 30
random.seed(20261111)
BIN = sys.argv[1]
decs = [(1e8, 1e9), (1e9, 1e10), (1e10, 1e11), (1e11, 5e11)]
ns = []
for lo, hi in decs:
    for _ in range(20):
        t = mpmath.mpf(lo) * (mpmath.mpf(hi) / lo) ** random.random()
        ns.append(int(mpmath.floor(mpmath.siegeltheta(t) / mpmath.pi)))
out = subprocess.run([BIN, '-check'] + [str(n) for n in ns], capture_output=True, text=True).stdout.split('\n')
rows = []
for line in out:
    if not line.startswith('CHECK'): continue
    d = dict(kv.split('=') for kv in line.split()[1:])
    rows.append(d)
print('decade  n  |dg| |dZ_new| |dZp_new| |dZ_old| |dZp_old| |dmu_new|/mu |dmu_old|/mu  (vs mpmath)')
stats = {}
for d in rows:
    n = int(d['n']); g = mpmath.grampoint(n)
    tnew = mpmath.mpf(d['tc']) + mpmath.mpf(d['d'])
    Zm = mpmath.siegelz(g); Zpm = mpmath.siegelz(g, derivative=1)
    Zn, Zpn, Zo, Zpo = (float(d[k]) for k in ('Z', 'Zp', 'old_Z', 'old_Zp'))
    mum = abs(Zpm / Zm)
    e = [abs(tnew - g), abs(Zn - Zm), abs(Zpn - Zpm), abs(Zo - Zm), abs(Zpo - Zpm), abs(abs(Zpn / Zn) - mum) / mum, abs(abs(Zpo / Zo) - mum) / mum]
    dec = int(mpmath.floor(mpmath.log10(g)))
    s = stats.setdefault(dec, [0.0] * 7 + [0])
    for i in range(7): s[i] = max(s[i], float(e[i]))
    s[7] += 1
    print(dec, n, ' '.join(f'{float(x):.1e}' for x in e), f'Z={float(Zm):.6f} Zp={float(Zpm):.4f}', flush=True)
print('MAX per decade (n points): |dg| |dZ_new| |dZp_new| |dZ_old| |dZp_old| rel_dmu_new rel_dmu_old')
for dec in sorted(stats):
    s = stats[dec]; print(f'1e{dec}: ({s[7]})', ' '.join(f'{x:.1e}' for x in s[:7]))
