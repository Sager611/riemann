import re,sys,collections
f = "/Users/adriansager/playground/riemann/.scratch/alp/lean/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean"
L = open(f, encoding="utf-8").read().split("\n")
decls = []
cur = None
for i, line in enumerate(L, 1):
    m = re.match(r"^(theorem|def|structure|lemma) (\S+)", line)
    if m:
        cur = [i, m.group(2), []]; decls.append(cur)
    elif cur is not None:
        cur[2].append(line)
kinds = collections.Counter()
uses_O_field = []
for i, name, body in decls:
    b = " ".join(x.strip() for x in body)
    if "OpenLiteratureStatements" in name: continue
    if re.search(r"\bO\.(caraianiLocalGlobalDual|localLanglandsExists|galoisRepOfPolarizedPairLit)", b) or re.search(r"\bO\.[0-9]", b):
        uses_O_field.append((i, name.split(".")[-1], b[:200]))
    head = b.split()[0] if b.split() else ""
    if re.search(r"\(O\)|\bO\b", b): kinds["mentions O"] += 1
    else: kinds["no O"] += 1
print("decls:", len(decls), dict(kinds))
print("direct field uses of O:")
for u in uses_O_field: print(" ", u[0], u[1], "::", u[2])
# sample of bodies that do not mention O
n=0
for i, name, body in decls:
    b = " ".join(x.strip() for x in body)
    if not re.search(r"\bO\b", b) and n < 6:
        print("noO", i, name.split(".")[-1], "::", b[:140]); n+=1
