
import re
mods = ["Arch", "Window", "Parity", "Det", "Num", "Exp", "Trig", "DetTrig"]
want = {"Arch": ["theoremE1P", "kappa2_gt_kstar", "e1p_enclosures", "inputs_consistentP", "nu_edge_iff", "nu_edge_eq", "exponent_edge"],
        "Window": ["theoremE1P_gauss", "H1P_of_gauss", "eta_edge", "eta_step", "gauss_consistent"],
        "Parity": ["theoremE1P_parity", "theoremE1P_fine", "H3P_of_inputs"],
        "Det": ["detector_threshold", "kstar_det_bounds", "kstar_anti", "cW_detp_ge", "cW_detm_le"],
        "Num": ["cW_edge_ge", "cW_edge_le", "kappa2_bounds", "kappa2_upper", "h_edge_gt", "h_edge_le", "kstar_edge_bounds"],
        "Exp": ["exp_Lhi2", "exp_Llo2", "sq_lo", "sq_up"], "Trig": ["base_sin", "cos_two", "trig_xlo2", "trig_xhi2"], "DetTrig": ["trig_det_lo", "trig_det_hi"]}
out = []
for m in mods:
    L = open(f"Riemann/Deep/B3Cert/Lean/{m}.lean").read().splitlines()
    for i, l in enumerate(L):
        mm = re.match(r"theorem (\S+)", l)
        if mm and mm.group(1) in want[m]:
            s = []
            j = i
            while True:
                t = L[j]; s.append(t)
                if t.rstrip().endswith(":= by") or t.rstrip().endswith(":=") or ":= by" in t or re.search(r":=\s*$", t) or (":=" in t and "⟨" in t):
                    break
                j += 1
            s[-1] = re.sub(r"\s*:=.*$", "", s[-1])
            out.extend(s)
defs = []
for m in ["Arch", "Window"]:
    L = open(f"Riemann/Deep/B3Cert/Lean/{m}.lean").read().splitlines()
    for i, l in enumerate(L):
        if re.match(r"(noncomputable )?def ", l):
            j = i; s = []
            while j < len(L) and L[j].strip() != "" and not L[j].startswith("/--"):
                s.append(L[j]); j += 1
            defs.extend(s)
open("Riemann/Deep/B3Cert/Lean/audit/stmts.txt", "w").write("\n".join(out) + "\n")
open("Riemann/Deep/B3Cert/Lean/audit/defs.txt", "w").write("\n".join(defs) + "\n")
print(len(out), len(defs))

