# Build the certification list: every CAND with mu < MU (screened) and every window/run minimum (WIN/RESO lines).
import sys, re
MU = float(sys.argv[1]); files = sys.argv[2:]
ns = {}
for fn in files:
    for l in open(fn):
        if l.startswith('CAND '):
            mu = float(re.search(r'mu=(\S+)', l).group(1)); n = int(re.search(r' n=(\d+)', l).group(1))
            if mu < MU: ns[n] = (mu, l.split()[1] + ' cand')
        elif l.startswith('WIN '):
            m = re.search(r'minmu=(\S+) n=(\d+)', l); ns.setdefault(int(m.group(2)), (float(m.group(1)), l.split()[1] + ' winmin'))
            m = re.search(r'minnorm=(\S+) n=(\d+) \(mu=(\S+)\)', l); ns.setdefault(int(m.group(2)), (float(m.group(3)), l.split()[1] + ' winminnorm'))
for n in sorted(ns, key=lambda k: ns[k][0]):
    print(n, ns[n][0], ns[n][1])
