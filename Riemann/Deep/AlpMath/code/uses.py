import sys, subprocess, collections
imp = collections.defaultdict(list)
for line in open("/Users/adriansager/playground/riemann/.scratch/alpmath/imports.txt", encoding="utf-8", errors="replace"):
    p, _, rest = line.partition(":import ")
    if not rest: continue
    imp[rest.strip()].append(p[2:] if p.startswith("./") else p)
base = "/Users/adriansager/playground/riemann/.scratch/alp/lean/src/"
for mod, name in [(a.split("@")[0], a.split("@")[1]) for a in sys.argv[1:]]:
    users = imp.get(mod, [])
    hits = []
    if users:
        r = subprocess.run(["rg", "-l", "-w", name] + [base + u for u in users], capture_output=True, text=True)
        hits = [h.replace(base, "") for h in r.stdout.split()]
    print(f"{mod.split('.')[-1]}::{name}: importers={len(users)} using={len(hits)} e.g. {hits[:3]}")
