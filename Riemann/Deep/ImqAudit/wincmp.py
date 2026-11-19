# wincmp.py M_FILE T4_FILE N: own exact M(d) (hcount MONLY, every fundamental d of the window) against isieve2's tier 4:
# {d : M(d) <= N} must equal the tier-4 set, with the same M values.
import sys
own = {}
for l in open(sys.argv[1]):
    a = l.split(); own[int(a[0])] = int(a[4])
N = int(sys.argv[3])
t4 = {int(l.split()[0]): int(l.split()[1]) for l in open(sys.argv[2])}
low = {d: m for d, m in own.items() if m <= N}
miss = [d for d in low if d not in t4]          # true M <= N but certified: tier 1 unsound
extra = [d for d in t4 if d not in low]         # tier 4 but true M > N (or not fundamental)
mdiff = [d for d in t4 if d in own and own[d] != t4[d]]
near = sorted(m - N for m in own.values() if m > N)[:3]
print('N=%d window fundamental d %d; own M<=N %d; tier4 %d; certified with M<=N (must be 0) %d; extra %d; M mismatch %d; smallest M-N above N %s'
      % (N, len(own), len(low), len(t4), len(miss), len(extra), len(mdiff), near))
