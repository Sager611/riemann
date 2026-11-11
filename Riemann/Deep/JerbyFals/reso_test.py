# Test the resonance proxy: rank of known low-mu candidates among proxy bins of a scanned window range.
import subprocess, sys, re
BIN = sys.argv[1]; log = sys.argv[2]; wid = sys.argv[3]; h = sys.argv[4]; SEP = float(sys.argv[5]); K = int(sys.argv[6])
w = [l for l in open(log) if l.startswith('WIN ' + wid + ' ')][0]
t0, t1 = map(float, re.search(r't=\[([^,]+),([^\]]+)\]', w).groups())
cands = [(float(re.search(r' t=(\S+)', l).group(1)), float(re.search(r'norm=(\S+)', l).group(1))) for l in open(log) if l.startswith('CAND ' + wid + ' ')]
p = subprocess.run([BIN, f'{t0:.6f}', f'{t1 - t0:.6f}', h, sys.argv[7], str(SEP), str(K), '40', 'x'], capture_output=True, text=True)
tops = [float(re.search(r' t=(\S+)', l).group(1)) for l in p.stderr.split('\n') if ' t=' in l]
Fs = [float(re.search(r' F=(\S+)', l).group(1)) for l in p.stderr.split('\n') if ' F=' in l]
nb = (t1 - t0) / SEP
print(f'{wid}: range {t1-t0:.0f}, bins {nb:.0f}, top {K} = {100*K/nb:.2f}% ; F range of top bins {Fs[-1]:.3f}..{Fs[0]:.3f}')
for tc, nm in sorted(cands, key=lambda x: x[1]):
    r = [k for k, tt in enumerate(tops) if abs(tt - tc) < 1.5 * SEP]
    print(f'  cand t={tc:.3f} norm={nm:.4f} -> proxy rank {r[0] if r else ">"+str(K)}')
