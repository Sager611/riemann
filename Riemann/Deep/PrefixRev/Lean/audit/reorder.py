# Reorder LEAN.md sections into the brief's order and append build/scan/axioms sections.
import re
p = "Riemann/Deep/PrefixRev/LEAN.md"
s = open(p).read()
head, rest = s.split("\n## ", 1)
secs = ("## " + rest).split("\n## ")
secs = [x if x.startswith("## ") else "## " + x for x in secs]
by = {}
for x in secs:
    title = x.split("\n", 1)[0][3:].strip()
    by[title] = x.rstrip() + "\n"
extra = open(".scratch/prefixlean/md_build.txt").read().rstrip() + "\n"
for x in ("\n" + extra).split("\n## ")[1:]:
    title = x.split("\n", 1)[0].strip()
    by[title] = "## " + x.rstrip() + "\n"
order = ["Files (Riemann/Deep/PrefixRev/Lean/, 17 files)",
 "Main Lean statements (verbatim; @-entries are Lean's elaborated #check output, which shows section hypotheses)",
 "Named hypotheses (Props passed as arguments; no axiom declarations)",
 "Build", "Static scan", "#print axioms", "Fidelity", "Memory cap and incidents", "Review", "Not done", "FORMALIZATIONS.md entry"]
missing = [t for t in order if t not in by]
extra_t = [t for t in by if t not in order]
out = head.rstrip() + "\n\n" + "\n".join(by[t] for t in order if t in by)
open(p, "w").write(out)
print("missing", missing, "unplaced", extra_t, "lines", out.count("\n"))
