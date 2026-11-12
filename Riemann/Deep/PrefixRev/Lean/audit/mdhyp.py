import re
blocks = open(".scratch/prefixlean/statements.txt").read().split("\n\n-- ")
byname = {}
for b in blocks:
    b = b if b.startswith("-- ") else "-- " + b
    m = re.search(r"^(theorem|def) (\S+)", b, re.M)
    if m: byname[m.group(2)] = "\n".join(l for l in b.strip().split("\n") if not l.startswith("-- "))
src = {
"OffLineLZZW": "Source: Lei, Zhao, Zhu, Wu, arXiv 2610.30892v1, Theorem 3.1, first sentence (TeX l. 221-223: \"If m+k<n-1, then <r_n,r_m,r_k> < S_n\"), standing assumption 2 <= k < m < n (TeX l. 109-111). Published, not formalized here (its graph/edge-count proof is NOTE section 1). Lean's natural n - 1 is exact (n >= 4).",
"Prop5Core": "Source: OUR Proposition 5, Steps 1-4 (NOTE section 4), restricted to the odd coprime line. OPEN OBLIGATION (new mathematics, not formalized; verify.py checks it for n <= 101). StarD is the NOTE's (★_d) in psi = phi - 1 coordinates.",
"JonesK0": "Source: G. A. Jones, Primitive permutation groups containing a cycle, Bull. Aust. Math. Soc. 89 (2014) 159-165, Theorem 1.2, case k = 0 (arXiv 1209.5169, TeX l. 54-72; rests on Jones, J. Group Theory 5 (2002) and CFSG). Stated with weakenings that Jones implies: (a) C_p <= G <= AGL_1(p) replaced by its consequence nfix <= 1 for nonidentity elements; (b) through the abstract predicate Coll q d H, intended meaning PGL_d(q) <= H <= PΓL_d(q) on the points of PG(d-1,q) (Mathlib: Projectivization; not instantiated: too heavy for the 290 MB cap); (c) G = L_2(11), M_11, M_23 replaced by the fixed-point bounds 3 (degree 11, M_11 sharply 4-transitive) and 7 (M_23 on 23 points).",
"Lemma7": "Source: NOTE Lemma 7 (fixed points of collineations; proof in NOTE section 5, brute-force check in verify.py). Named, not proved here.",
"Prop5Primitive": "Used by Main.theoremA (older interface): primitivity on the odd line. Proved from Prop5Core by PrimB2.prop5_of_core (Lemma 3 + (★_d)).",
"JonesLemma7": "Used by Main.theoremA (older interface): the merged consequence of JonesK0 + Lemma7. theoremA_jones (Jones.lean) uses JonesK0 and Lemma7 separately instead; prefer it.",
}
out = ["## Named hypotheses (Props passed as arguments; no axiom declarations)", ""]
for nm in ["OffLineLZZW", "Prop5Core", "JonesK0", "Lemma7", "Prop5Primitive", "JonesLemma7"]:
    out += ["### " + nm, "", src[nm], "", "~~~lean", byname[nm], "~~~", ""]
out += ["Main entry point: theoremA_jones (hOff : OffLineLZZW) (hC : Prop5Core) (hJ : JonesK0 Coll) (hL : Lemma7 Coll). theoremA_necessary needs only OffLineLZZW; line_necessary needs nothing.", ""]
open(".scratch/prefixlean/md_hyp.txt", "w").write("\n".join(out))
print(len(out))
