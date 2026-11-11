# ctable.py: counting function N(kappa) from laurent_count.py outputs (C_*.json), by S and tmax.
import sys, json, glob
KS = ['0.5', '0.3', '0.2', '0.15', '0.1', '0.07', '0.05', '0.04', '0.03', '0.02']
print('%-22s %3s %6s %5s |' % ('file', '|S|', 'tmax', 'N') + ' '.join('%6s' % k for k in KS))
for f in sys.argv[1:]:
    d = json.load(open(f)); I = d['info']
    print('%-22s %3d %6g %5d |' % (f.split('/')[-1][:-5], len(I['primes']) + 1, I['tmax'], I['N']) + ' '.join('%6d' % d['N'][k] for k in KS))
