
import re
V = [l for l in open("Riemann/Deep/B3Cert/Lean/audit/verify.log").read().splitlines() if " rc=" in l]
P = [l for l in open("Riemann/Deep/B3Cert/Lean/audit/peaks.log").read().splitlines() if " rc=" in l]
rows = ["| module | wall | footprint (time -l) | memcap peak (top MEM) |", "|---|---|---|---|"]
for l in V:
    m = re.match(r"Riemann/Deep/(\w+)/Lean/(\w+) rc=(\d+) (\d+)s footprint=(\d+)MiB memcap: peak (\d+) MB", l)
    d, mod, rc, t, fp, mc = m.groups()
    tag = "E1Lean " if d == "SantiAudit" else ""
    rows.append(f"| {tag}{mod} | {t} s | {fp} MiB | {mc} MB{' (exited between samples)' if mc == '0' else ''} |")
allfp = [int(x) for x in re.findall(r"footprint=(\d+)MiB", "\n".join(V + P))]
allmc = [int(x) for x in re.findall(r"memcap: peak (\d+) MB", "\n".join(V + P))]
own = [l for l in V + P if "B3Cert" in l]
ownfp = [int(x) for x in re.findall(r"footprint=(\d+)MiB", "\n".join(own))]
ownmc = [int(x) for x in re.findall(r"memcap: peak (\d+) MB", "\n".join(own))]
out = "\n".join(rows)
out += f"\n\nAll compiles of the session ({len(V) + len(P)} runs, peaks.log + verify.log): footprint ≤ {max(allfp)} MiB, top MEM ≤ {max(allmc)} MB (both E1Lean's Parity); E1PLean files alone: footprint ≤ {max(ownfp)} MiB, top MEM ≤ {max(ownmc)} MB; one failed compile (Num, a reduced-literal mismatch with Exp's c_hi, 125 s, 262 MB; fixed by writing c_hi unreduced). Scripts (Python, exact Fractions): ≤ 1 s each, below memcap's first sample. Every run under memcap 290 MB; all below the 300 MB cap."
open("Riemann/Deep/B3Cert/Lean/audit/build_table.md", "w").write(out + "\n")
print(out[-700:])

