import os, re, sys, json, collections
# collect declaration full names, namespace-aware (approximate), from .lean files under roots
roots = sys.argv[1].split(',')
targets = sys.argv[2].split(',')   # short names to look for, e.g. dedekindZeta,DifferentiableOn
prefix = sys.argv[3] if len(sys.argv) > 3 else ''
decl_re = re.compile(r'^\s*(?:@\[[^\]]*\]\s*)*(?:(?:private|protected|noncomputable|public|nonrec|partial|unsafe|scoped|local)\s+)*(def|theorem|lemma|abbrev|instance|structure|class|inductive|opaque|axiom|irreducible_def)\s+([^\s:({\[]+)')
ns_re = re.compile(r'^\s*namespace\s+(\S+)')
sec_re = re.compile(r'^\s*(?:noncomputable\s+|public\s+)*section\b\s*(\S*)')
end_re = re.compile(r'^\s*end\b\s*(\S*)')
strip_comments = re.compile(r'/-.*?-/', re.S)
hits = []; ndecl = 0
for root in roots:
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d not in ('.lake', '.git')]
        for f in fn:
            if not f.endswith('.lean'): continue
            p = os.path.join(dp, f)
            src = open(p, encoding='utf-8', errors='replace').read()
            src = strip_comments.sub(lambda m: '\n' * m.group(0).count('\n'), src)
            stack = []  # entries: ('ns', name) or ('sec', name)
            for ln, line in enumerate(src.split('\n'), 1):
                line = line.split('--')[0]
                m = ns_re.match(line)
                if m: stack.append(('ns', m.group(1))); continue
                m = sec_re.match(line)
                if m and not line.strip().startswith('sectionX'): stack.append(('sec', m.group(1))); continue
                m = end_re.match(line)
                if m:
                    if stack: stack.pop()
                    continue
                m = decl_re.match(line)
                if m:
                    ndecl += 1
                    name = m.group(2)
                    nss = '.'.join(n for k, n in stack if k == 'ns')
                    full = name[len('_root_.'):] if name.startswith('_root_.') else ((nss + '.' + name) if nss else name)
                    last = full.split('.')[-1]
                    for t in targets:
                        tparts = t.split('.')
                        if full.endswith('.' + t) or full == t or (len(tparts) == 1 and last == t):
                            hits.append((full, m.group(1), p, ln))
print('declarations seen:', ndecl)
for full, kind, p, ln in sorted(set(hits)):
    if prefix and not full.startswith(prefix): continue
    print(kind, full, p.split('/hpaudit/')[-1], ln)
