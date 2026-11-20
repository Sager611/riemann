import sys, collections
src = "/Users/adriansager/playground/riemann/.scratch/alp/lean/src/"
ids = {}; names = []
def gid(m):
    i = ids.get(m)
    if i is None:
        i = len(names); ids[m] = i; names.append(m); adj.append([])
    return i
adj = []
with open("/Users/adriansager/playground/riemann/.scratch/alpmath/imports.txt", encoding="utf-8", errors="replace") as f:
    for line in f:
        p, _, rest = line.partition(":import ")
        if not rest: continue
        if p.startswith("./"): p = p[2:]
        if not p.endswith(".lean"): continue
        mod = p[:-5].replace("/", ".")
        tgt = rest.strip().split()[0] if rest.strip() else None
        if tgt is None: continue
        a = gid(mod); b = gid(tgt)
        adj[a].append(b)
def cone(root):
    r = ids[root]; seen = bytearray(len(names)); seen[r] = 1; st = [r]
    while st:
        x = st.pop()
        for y in adj[x]:
            if not seen[y]:
                seen[y] = 1; st.append(y)
    return seen
for root in ["Main", "Supplements.Kernels.MainTheorems", "Development.AutomorphicArithmetic.CharacterBoxes.MainConclusion", "Formalism.GammaFactorNonvanishing.CritLine", "Development.AutomorphicArithmetic.CharacterBoxes.RiemannHypothesisOfEveryCharacter"]:
    s = cone(root)
    tot = sum(s); local = sum(1 for i in range(len(names)) if s[i] and not names[i].startswith("Mathlib") and names[i] != "Mathlib")
    tops = collections.Counter(names[i].split(".")[0] for i in range(len(names)) if s[i])
    print(root, "cone files:", tot, "non-Mathlib:", local, dict(tops.most_common(10)))
    if root == "Supplements.Kernels.MainTheorems":
        with open("/Users/adriansager/playground/riemann/.scratch/alpmath/cone_MT.txt", "w") as g:
            for i in range(len(names)):
                if s[i] and not names[i].startswith("Mathlib"): g.write(names[i] + "\n")
print("modules seen:", len(names))
# SmallClassNumbers in MT cone?
s = cone("Supplements.Kernels.MainTheorems")
for m in ["Development.AutomorphicArithmetic.AutomorphicVarieties.SmallClassNumbers.CorGKOf", "Development.AutomorphicArithmetic.AutomorphicVarieties.SmallClassNumbers.ClassNumberThree", "Development.AutomorphicArithmetic.AutomorphicVarieties.SmallClassNumbers.OnePoint", "Development.AutomorphicArithmetic.AutomorphicVarieties.SmallClassNumbers.LenEqZeroIff"]:
    print(m.split(".")[-1], "in MT cone:", bool(s[ids[m]]) if m in ids else "unknown")
