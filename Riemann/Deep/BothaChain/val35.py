# val35.py (bothachain): validation of the precision change (bits factor 1.6 -> 1.25, mkint.py BITSF) at the overlap point c = 35 against
# CCMFMM's certified window (logs/int.jsonl: U(35) = [3.35474138011541e-177 +/- 1.10e-192] at NG = 320, 1024 bits; odd certified at level
# 8.972773e-175, nu = 8.972773e-174, kappa eta = 3.566e-269). Here: gstate at NG = 320 with BITS bits (U must agree in all printed digits) and the
# odd certificate at the same level with nu = level/10 (our chain's rung) at BITS bits (must pass; kappa eta printed).   usage: val35.py <bits>
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, 'ccm'))
from run_ccm import gstate, cert
bits = int(sys.argv[1]); t0 = time.time()
g, gm = gstate('35', 320, bits, 0, 'g35_%d' % bits)
ro = cert('35', 'odd', '8.972773e-175', bits, 2043, ['8.972773e-176'], 'v35o_%d' % bits)
print(json.dumps(dict(c='35', bits=bits, U=gm['U'], U_ccmfmm='[3.35474138011541e-177 +/- 1.10e-192]', g_sec=gm['sec'], g_rss=gm['rss_mb'],
                      odd=dict((k, ro.get(k)) for k in ('rigorous', 'nu_cert', 'kappa_eta', 'n', 'nU', 'rss_mb', 'sec')), sec=round(time.time() - t0, 1))))
