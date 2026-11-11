
import os, re, json, sys
ROOT = "/Users/adriansager/playground/riemann"
def path_of(mod):
    return os.path.join(ROOT, *mod.split(".")) + ".lean"
imp_re = re.compile(r'^\s*(?:public\s+|meta\s+|private\s+)*import\s+(\S+)', re.M)
cache = {}
def imports(mod):
    if mod in cache: return cache[mod]
    p = path_of(mod)
    if not os.path.exists(p):
        cache[mod] = None; return None
    src = open(p, encoding="utf-8").read()
    # header only: imports before first non-import, non-comment line
    res = []
    for m in imp_re.finditer(src):
        res.append(m.group(1))
    cache[mod] = res
    return res
target = "Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Energy.CertifiedExistence"
closure = set(); stack=[target]; external=set()
while stack:
    m = stack.pop()
    if m in closure: continue
    if not m.startswith("Riemann."):
        external.add(m); continue
    imps = imports(m)
    if imps is None:
        print("MISSING", m); continue
    closure.add(m)
    stack.extend(imps)
print("closure (Riemann.*):", len(closure))
print("external roots:", len(external), sorted(external)[:20])
pg = "Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Moments.PrimeGlobal"
print("PrimeGlobal in closure:", pg in closure)
# downstream of PrimeGlobal within closure
down = {pg}
changed=True
while changed:
    changed=False
    for m in closure:
        if m in down: continue
        if any(i in down for i in imports(m)):
            down.add(m); changed=True
print("downstream incl PG:", len(down))
lines = sum(open(path_of(m),encoding='utf-8').read().count("\n") for m in down)
print("lines:", lines)
from collections import Counter
c = Counter(".".join(m.split(".")[5:6]) for m in down)
print(c.most_common(10))
# check: no file in closure outside down imports a file in down (by def of down, trivially true); check entire repo
allmods=[]
for dp,dn,fn in os.walk(os.path.join(ROOT,"Riemann")):
    for f in fn:
        if f.endswith(".lean"):
            allmods.append(os.path.relpath(os.path.join(dp,f),ROOT)[:-5].replace("/","."))
outside_importers=[m for m in allmods if m not in down and any(i in down for i in (imports(m) or []))]
print("repo files outside down importing down:", len(outside_importers), outside_importers[:10])
up = closure - down
print("upstream (to build):", len(up))
frontier = sorted({i for m in down for i in imports(m) if i in up})
print("frontier size:", len(frontier))
# frontier maximal elements: those in up not imported by any other in up? we just need set of imports from down into up
json.dump({"closure":sorted(closure),"down":sorted(down),"up":sorted(up),"frontier":frontier,"external":sorted(external)}, open("/Users/adriansager/playground/riemann/.scratch/oaisext/closure.json","w"), indent=0)
# minimal set of targets covering up: elements of up not imported by any other element of up
imported_by_up = {i for m in up for i in imports(m)}
tops = sorted(m for m in up if m not in imported_by_up)
print("up tops:", len(tops))
for t in tops: print("  ", t)

