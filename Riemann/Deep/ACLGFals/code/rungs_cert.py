# rungs_cert.py: exact delta_4, delta_6, delta_8 at the windows certified by Chekmachev (Zenodo 23848597, audited in ../../ChekAudit)
# and weilsimple, for the certified LA/gap ratios. Prints JSON.
import sys, os, json, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prolate import rungs
out = {}
for c in [int(x) for x in sys.argv[1:]]:
    R, dps, K = rungs(mp.mpf(c), 10)
    out[c] = {n: mp.nstr(R[n]['delta'], 12) for n in (4, 6, 8, 10)}
print(json.dumps(out))
