import re
rows = []
for l in open(".scratch/prefixlean/peaks.log"):
    p = l.split()
    if len(p) < 3: continue
    m = p[1]
    fp = re.search(r"(\d+) peak memory footprint", l)
    mc = re.search(r"memcap: peak (\d+) MB", l)
    rt = re.search(r"([\d.]+) real", l)
    rc = re.search(r" rc=(\d+)", l)
    rows.append((p[0], m, rc.group(1), int(fp.group(1)) / 1e6 if fp else 0, mc.group(1) if mc else "-", rt.group(1) if rt else "-"))
out = ["## Build", "",
"Clean rebuild 2026-11-12 18:18-18:57 UTC (.scratch/prefixlean/verify.sh; Main's 1650 s wall time is a host stall, 4.5 s user). One module at a time, no lake process, oleans to .scratch/prefixlean/olean (build output, deleted at the end):", "",
"    cp Riemann/Deep/DPCUncond/Lean/audit/leanpath.txt .scratch/prefixlean/leanpath.txt   # lake's LEAN_PATH",
"    # for m in Defs CycleA CycleB CycleC Parity IsCycleP Small BlocksA BlocksB Main PrimA PrimB1 PrimB2 Jones Orbits CheckA CheckB:",
"    LEAN_PATH=.scratch/prefixlean/olean:$(cat .scratch/prefixlean/leanpath.txt) bash .scratch/memcap.sh 290 1800 \\",
"      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 \\",
"      -o .scratch/prefixlean/olean/Riemann/Deep/PrefixRev/Lean/$m.olean Riemann/Deep/PrefixRev/Lean/$m.lean", "",
"| Module | rc | peak footprint (time -l) | memcap peak | wall |", "|---|---|---|---|---|"]
for t, m, rc, fp, mc, rt in rows:
    out.append("| %s | %s | %.0f MB | %s MB | %s s |" % (m, rc, fp, mc, rt))
out += ["", "All 17 modules rc = 0; maximum footprint %.0f MB (Small). Kernel checks: decide +kernel on Bool evaluations (Small.lean)." % max(r[3] for r in rows), ""]
out += ["## Static scan", "", "    rg -n -w 'sorry|admit|native_decide|skipKernelTC' Riemann/Deep/PrefixRev/Lean   # no matches",
        "    rg -n '^\\s*axiom\\s' Riemann/Deep/PrefixRev/Lean                              # no matches", ""]
ax = []
for f in [".scratch/prefixlean/logs/CheckA.log", ".scratch/prefixlean/logs/CheckB.log"]:
    for l in open(f):
        if "depends on axioms" in l: ax.append(l.rstrip())
out += ["## #print axioms", "", "CheckA.lean (40) and CheckB.lean (9), clean rebuild:", "", "~~~"] + ax + ["~~~", ""]
open(".scratch/prefixlean/md_build.txt", "w").write("\n".join(out))
print(len(ax), len(rows))
