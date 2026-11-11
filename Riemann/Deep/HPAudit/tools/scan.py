#!/usr/bin/env python3
"""Red-flag scanner for Lean sources. Usage: scan.py OUTPREFIX ROOT [ROOT...]
Strips comments (nested block comments, line comments) and string literals, then
reports pattern hits. Writes OUTPREFIX.hits.txt and prints a summary."""
import os, re, sys, collections

PATTERNS = collections.OrderedDict([
  ("sorry", r"\bsorry\b"),
  ("admit", r"\badmit\b"),
  ("axiom_decl", r"(^|\s)(private\s+|protected\s+|noncomputable\s+)*axiom\s"),
  ("axioms_kw", r"\baxioms?\b"),
  ("native_decide", r"\bnative_decide\b"),
  ("decide_kernel_cfg", r"decide\s*\+\s*kernel|kernel\s*:=\s*true|decide\s*:=\s*true"),
  ("ofReduceBool", r"ofReduceBool|reduceBool|ofReduceNat|reduceNat"),
  ("implemented_by", r"implemented_by"),
  ("extern", r"\bextern\b"),
  ("unsafe", r"\bunsafe\b"),
  ("opaque", r"\bopaque\b"),
  ("partial", r"\bpartial\s+def\b|\bpartial\b"),
  ("skipKernelTC", r"skipKernelTC|debug\.skip"),
  ("addDecl", r"addDeclWithoutChecking|addDecl\b|Environment\.addDecl|replaceConst|modifyEnv|setEnv"),
  ("Environment", r"\bEnvironment\b"),
  ("run_cmd", r"\brun_cmd\b|\brun_elab\b|\brun_meta\b"),
  ("elab", r"(^|\s)elab\b|\belab_rules\b|\bmacro_rules\b|(^|\s)macro\b|(^|\s)syntax\b"),
  ("initialize", r"\binitialize\b|\bbuiltin_initialize\b"),
  ("hash_eval", r"#eval|#exit|#guard_msgs|#reduce"),
  ("set_option_kernel", r"set_option\s+\S*kernel\S*|set_option\s+debug\S*|set_option\s+trust\S*|set_option\s+compiler\S*"),
  ("trustCompiler", r"trustCompiler|Lean\.trust"),
  ("csimp", r"@\[\s*csimp|csimp\b"),
  ("import_Lean", r"^\s*import\s+Lean\b"),
  ("notation", r"\bnotation\b|\binfix[lr]?\b|\bprefix\b|\bpostfix\b"),
  ("instance_prio", r"instance\s*\(\s*priority"),
  ("local_riemannZeta_def", r"\b(def|abbrev|noncomputable def)\s+(riemannZeta|LFunction|dedekindZeta)\b"),
  ("print_axioms", r"#print\s+axioms"),
])
COMPILED = {k: re.compile(v, re.M) for k, v in PATTERNS.items()}

def strip(src):
    out = []; i = 0; n = len(src); depth = 0
    while i < n:
        if depth == 0 and src.startswith("--", i):
            j = src.find("\n", i); j = n if j < 0 else j
            i = j; continue
        if src.startswith("/-", i):
            depth += 1; i += 2; continue
        if depth > 0:
            if src.startswith("-/", i): depth -= 1; i += 2; continue
            if src[i] == "\n": out.append("\n")
            i += 1; continue
        c = src[i]
        if c == '"':
            j = i + 1
            while j < n and src[j] != '"':
                if src[j] == "\\": j += 1
                j += 1
            out.append('""'); i = j + 1; continue
        if c == "'" and i + 2 < n and src[i+2] == "'" and src[i+1] != "\\":
            out.append("' '"); i += 3; continue
        out.append(c); i += 1
    return "".join(out)

def main():
    prefix = sys.argv[1]; roots = sys.argv[2:]
    hits = collections.defaultdict(list); nfiles = 0; nlines = 0
    for root in roots:
        for dp, dn, fn in os.walk(root):
            dn[:] = [d for d in dn if d not in (".lake", ".git")]
            for f in fn:
                if not f.endswith(".lean"): continue
                p = os.path.join(dp, f)
                src = open(p, encoding="utf-8", errors="replace").read()
                nfiles += 1; nlines += src.count("\n") + 1
                s = strip(src)
                lines = s.split("\n")
                for k, rx in COMPILED.items():
                    for m in rx.finditer(s):
                        ln = s.count("\n", 0, m.start()) + 1
                        hits[k].append((p, ln, lines[ln-1].strip()[:200]))
    with open(prefix + ".hits.txt", "w") as fh:
        fh.write(f"files {nfiles} lines {nlines}\n")
        for k in PATTERNS:
            fh.write(f"== {k}: {len(hits[k])}\n")
            for p, ln, t in hits[k]:
                fh.write(f"{p}:{ln}: {t}\n")
    print(f"files {nfiles} lines {nlines}")
    for k in PATTERNS:
        print(f"{k}: {len(hits[k])}")

main()

