import re
blocks = open(".scratch/prefixlean/statements.txt").read().split("\n\n-- ")
byname = {}
for b in blocks:
    b = b if b.startswith("-- ") else "-- " + b
    m = re.search(r"^(theorem|def) (\S+)", b, re.M)
    if m: byname[m.group(2)] = b.strip()
sig = open(".scratch/prefixlean/signatures.txt").read()
def sigof(name):
    m = re.search(r"(@PrefixRev\." + name + r" :.*?)(?=\n@PrefixRev\.|\ndef PrefixRev\.|\nmemcap)", sig, re.S)
    return m.group(1).strip() if m else None
order = [("Definitions (Defs, CycleA, Main, PrimB1, Jones, Small, Orbits)", ["revFun","rev","G","P","phi","IsBlockOf","IsPrimitive","StarD","nfix","Conjecture4Literal","cls"]),
 ("Lemma 1, Lemma 2, B.1 n-cycle", ["P_val","P_firstReturn","@phi_P","@phi_rn","@phi_rk","P_cycle","P_isCycle","P_support","G_eq_top_of_swap"]),
 ("Parity", ["sign_rev","G_le_alt_iff","eq_top_of_alt"]),
 ("Necessity on the line (B.1 intransitivity, B.2, B.3 blocks, parity)", ["G_ne_top_of_gcd","G_le_intervals","G_ne_top_of_dvd","G_ne_top_of_even","G_le_alt_of_mod4","line_necessary","theoremA_necessary","@G_le_cls","cls_surj"]),
 ("Counterexamples (kernel-checked)", ["G_7_4_2","G_11_6_4","G_13_10_2","not_conjecture4Literal"]),
 ("Lemma 3 and (★_d)", ["@block_classes","@star_of_block","prop5_of_core"]),
 ("Theorem A and B.4", ["theoremA_jones","line_B4_jones","theoremA","line_B4","theoremA_core"])]
out = ["## Main Lean statements (verbatim; @-entries are Lean's elaborated #check output, which shows section hypotheses)", ""]
for title, names in order:
    out.append("### " + title); out.append(""); out.append("~~~lean")
    for nm in names:
        if nm.startswith("@"):
            s = sigof(nm[1:]); out.append(s if s else "-- missing " + nm)
        else:
            b = byname.get(nm)
            if b is None: out.append("-- missing " + nm); continue
            lines = [l for l in b.split("\n") if not l.startswith("-- ")]
            out.append("\n".join(lines))
        out.append("")
    out.append("~~~"); out.append("")
open(".scratch/prefixlean/md_statements.txt", "w").write("\n".join(out))
print(len(out), sum(1 for l in out if "missing" in l))
