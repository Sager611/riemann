#!/usr/bin/env python3
# Generate TABLE.md section (a): one row per syntax/macro/elab/unif_hint declaration in the archive.
import re, sys
D, OUT = sys.argv[1], sys.argv[2]
rows = [l.rstrip('\n').split('\t') for l in open(D + 'syn91.tsv')]
special = {
 'kexactB': "tactic elab (Imports.lean:5-8): elabTermAndSynthesize e none, then g.assign without isDefEq; the elaborator's type check is skipped, the kernel still checks the finished declaration (an ill-typed assignment fails at addDecl)",
 'kcopy': "tactic elab (OwnStruct.lean:53-58): rebuilds own-structure values field by field (Expr.proj, mkConst ctor) and assigns unchecked; kernel checks the result",
 '#kcopy_probe': "command elab (OwnStruct.lean:60-64): diagnostic only, runs copyAt in probe mode (Kernel.isDefEq read-only) and logInfo; adds nothing",
 'cpt_copy_thms': "command elab (CutCptCopy.lean:68-91): re-adds a Formalism.* dependency cone with constants renamed (suffix Cpt / special map) via Lean.addDecl (.defnDecl/.thmDecl): kernel-checked (no debug.skipKernelTC anywhere)",
}
out = ['| # | file:line | scope | kind | atom | what it does | skips kernel? | changes a Main statement? | captures an identifier? |', '|---|---|---|---|---|---|---|---|---|']
for i, (sc, loc, txt) in enumerate(rows, 1):
    m = re.search(r'"([^"]*)"', txt); atom = (m.group(1).strip() if m else 'unif_hint')
    cat = re.search(r':\s*(tactic|term|command)\s*=>', txt)
    kind = 'elab ' if re.search(r'\belab\b', txt.split('"')[0]) else ('unif_hint ' if 'unif_hint' in txt else 'macro ')
    kind += (cat.group(1) if cat else ('-' if 'unif_hint' in txt else 'term'))
    if atom in special: what = special[atom]
    elif 'unif_hint' in txt: what = 'local unification hint for the elaborator; the kernel never sees hints'
    elif 'kexactB' in txt and 'term' in kind: what = 'term macro → (by kexactB e) (or a structure literal of such): unchecked goal assignment, kernel-checked declaration'
    elif 'tactic' in kind: what = 'tactic macro: expands to a tactic script; output is an ordinary proof term checked by the kernel'
    else: what = 'term macro: syntactic abbreviation, expanded before elaboration'
    reach = {'local': 'file-local', 'scoped': 'active only after open of its namespace', 'GLOBAL': 'global'}[sc]
    cap = 'no: new atom "' + atom + '" is none of Main.lean\'s tokens; ' + ('tactic category, Main.lean has no tactic block' if 'tactic' in kind else reach)
    stmt = 'no: ' + reach + '; Main.lean has no open, no namespace, no tactic' 
    out.append('| %d | %s | %s | %s | %s | %s | no | %s | %s |' % (i, loc, sc.lower(), kind, atom.replace('|', '\\|'), what, stmt, cap))
open(OUT, 'a').write('\n'.join(out) + '\n')
print(len(rows), 'rows')
