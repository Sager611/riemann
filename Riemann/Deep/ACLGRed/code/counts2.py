import json, math
rows = [json.loads(l) for l in open('/Users/adriansager/playground/riemann/.scratch/aclgred/spec_all.jsonl')]
sig = [2, 3, 4]
print("| λ² | N | dim/sector | ev<1e-10 | odd<1e-10 | ev<1e-2 | odd<1e-2 | " + " | ".join("#(1e-2, ½log%d)" % s for s in sig) + " | #(1e-2, 1) | first non-min ev / odd |")
print("|" + "---|"*(10 + len(sig) - 2))
for r in sorted(rows, key=lambda r: float(r['c'])):
    ev = sorted(float(x) for x in r['even']); od = sorted(float(x) for x in r['odd'])
    a = ev + od
    def cnt(lo, hi): return sum(1 for x in a if lo <= x < hi)
    nmev = [x for x in ev if x >= 1e-2]; nmod = [x for x in od if x >= 1e-2]
    cells = [r['c'], str(r['N']), "%d/%d" % (len(ev), len(od)), str(sum(1 for x in ev if x < 1e-10)), str(sum(1 for x in od if x < 1e-10)),
             str(sum(1 for x in ev if x < 1e-2)), str(sum(1 for x in od if x < 1e-2))]
    cells += [str(cnt(1e-2, 0.5*math.log(s))) for s in sig] + [str(cnt(1e-2, 1.0)), "%.3f / %.3f" % (nmev[0], nmod[0])]
    print("| " + " | ".join(cells) + " |")
    neg = [x for x in a if x < 0]
    if neg: print("NEGATIVE", r['c'], neg[:5])
