import re, glob, os
want = ["revFun","rev","G","P","phi","IsBlockOf","IsPrimitive","OffLineLZZW","Prop5Primitive","JonesLemma7","StarD","Prop5Core","nfix","JonesK0","Lemma7","Conjecture4Literal","cls",
 "P_val","P_firstReturn","phi_P","phi_rn","phi_rk","P_cycle","P_isCycle","P_support","G_eq_top_of_swap","sign_rev","G_le_alt_iff","eq_top_of_alt",
 "G_intransitive","G_ne_top_of_gcd","G_le_intervals","G_ne_top_of_dvd","G_ne_top_of_even","G_le_alt_of_mod4","line_necessary","theoremA_necessary",
 "line_alt_le","line_B4","theoremA","block_classes","star_of_block","prop5_of_core","theoremA_core","line_alt_le_jones","theoremA_jones","line_B4_jones",
 "G_le_cls","cls_surj","G_7_4_2","G_11_6_4","G_13_10_2","not_conjecture4Literal"]
out = []
d = "Riemann/Deep/PrefixRev/Lean"
for f in sorted(glob.glob(d + "/*.lean")):
    if "Check" in f: continue
    lines = open(f).read().split("\n")
    i = 0
    while i < len(lines):
        m = re.match(r"^(theorem|def) (\S+)", lines[i])
        if m and m.group(2) in want:
            # include preceding docstring
            j = i - 1
            doc = []
            while j >= 0 and (lines[j].startswith("/--") or lines[j].startswith("  ") or lines[j].endswith("-/")) and not lines[j].startswith("theorem") and lines[j].strip() != "":
                doc.insert(0, lines[j]); 
                if lines[j].startswith("/--"): break
                j -= 1
            if not (doc and doc[0].startswith("/--")): doc = []
            block = doc + [lines[i]]
            k = i + 1
            if m.group(1) == "def":
                while k < len(lines) and lines[k].strip() != "" and not lines[k].startswith("theorem") and not lines[k].startswith("/--"):
                    block.append(lines[k]); k += 1
            else:
                cur = lines[i]
                while ":= by" not in cur and not cur.rstrip().endswith(":=") and " := " not in cur and k < len(lines):
                    cur = lines[k]; block.append(cur); k += 1
            out.append("-- " + os.path.basename(f) + "\n" + "\n".join(block))
        i += 1
open(".scratch/prefixlean/statements.txt", "w").write("\n\n".join(out) + "\n")
print(len(out), sum(len(o) for o in out))
