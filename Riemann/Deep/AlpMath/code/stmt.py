import sys,re
base = "/Users/adriansager/playground/riemann/.scratch/alp/lean/src/"
for arg in sys.argv[1:]:
    f, ln = arg.rsplit(":", 1); ln = int(ln)
    L = open(base + f, encoding="utf-8").read().split("\n")
    out = []
    for i in range(ln - 1, min(ln + 14, len(L))):
        out.append(L[i])
        if ":=" in L[i]: break
    s = " ".join(x.strip() for x in out)
    s = s.split(":=")[0]
    print("@" + f.split("/")[-1] + ":" + str(ln) + " :: " + s[:700])
