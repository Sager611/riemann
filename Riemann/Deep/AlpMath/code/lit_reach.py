import re, collections
f = "/Users/adriansager/playground/riemann/.scratch/alp/lean/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean"
L = open(f, encoding="utf-8").read().split("\n")
decls = collections.OrderedDict(); cur = None; types = {}
for i, line in enumerate(L, 1):
    m = re.match(r"^(theorem|def|structure|lemma) (\S+)", line)
    if m:
        name = m.group(2).split(".")[-1]; cur = name; decls[cur] = [i, line]
        t = re.search(r"OpenLiteratureStatements\) : ([^ ]+)", line)
        types[name] = t.group(1).split(".")[-1] if t else ""
    elif cur: decls[cur].append(line)
body = {n: " ".join(v[1:]) for n, v in decls.items()}
root = "every_character_of_literature_of_vWin"
seen = {root}; st = [root]
while st:
    x = st.pop()
    for y in decls:
        if y not in seen and re.search(r"LiteratureTheorems\." + re.escape(y) + r"(?![A-Za-z0-9_'])", body[x]):
            seen.add(y); st.append(y)
print("decls", len(decls), "reachable from root", len(seen))
unreach = [n for n in decls if n not in seen and n != "OpenLiteratureStatements"]
print("unreachable:", len(unreach), [ (n, types.get(n)) for n in unreach][:40])
for t in ["GreenRegularCharacters", "HilbertClassFieldJ", "ComplexMultiplicationJSt", "RadfordShufflePolynomial", "SiegelKlingenDeligneRibet", "DeligneRibet", "UnusedSlot72", "Linnik", "WeylBranching", "QuadraticClassFieldTheory", "GreenCuspidal"]:
    ns = [n for n in decls if types.get(n) == t]
    print(t, "->", [(n, decls[n][0], n in seen) for n in ns])
