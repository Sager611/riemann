import os, re, sys, collections
roots = sys.argv[1].split(',')
closure = None
if len(sys.argv) > 2 and sys.argv[2]:
    import json; closure = set(json.load(open(sys.argv[2])))
strip_comments = re.compile(r'/-.*?-/', re.S)
inst_re = re.compile(r'^(\s*)(?:@\[[^\]]*\]\s*)*((?:(?:private|protected|noncomputable|public|scoped|local)\s+)*)instance\b(.*)$')
res = []
for root in roots:
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d not in ('.lake', '.git')]
        for f in fn:
            if not f.endswith('.lean'): continue
            p = os.path.join(dp, f)
            src = open(p, encoding='utf-8', errors='replace').read()
            src = strip_comments.sub(lambda m: '\n' * m.group(0).count('\n'), src)
            lines = src.split('\n')
            for i, line in enumerate(lines):
                m = inst_re.match(line.split('--')[0])
                if not m: continue
                mods = m.group(2)
                if 'local' in mods or 'scoped' in mods: continue
                # gather header until := or where
                hdr = m.group(3); j = i
                while not re.search(r':=|\bwhere\b', hdr) and j + 1 < len(lines) and j - i < 12:
                    j += 1; hdr += ' ' + lines[j].split('--')[0].strip()
                hdr = re.split(r':=|\bwhere\b', hdr)[0]
                # conclusion: text after the last top-level ':'
                depth = 0; last = -1
                for k, ch in enumerate(hdr):
                    if ch in '([{⦃': depth += 1
                    elif ch in ')]}⦄': depth -= 1
                    elif ch == ':' and depth == 0: last = k
                concl = hdr[last+1:].strip() if last >= 0 else hdr.strip()
                res.append((p, i + 1, concl[:160], hdr.strip()[:200]))
print('global instances:', len(res))
basic = re.compile(r'ℂ|ℝ|\bComplex\b|\bReal\b')
cls = collections.Counter()
for p, ln, concl, hdr in res:
    head = concl.split()[0] if concl else ''
    cls[head] += 1
print('instance classes (top 40):')
for k, v in cls.most_common(40): print('  ', v, k)
print('--- instances whose conclusion mentions ℂ/ℝ/Complex/Real:')
for p, ln, concl, hdr in res:
    if basic.search(concl):
        print(p.split('/hpaudit/')[-1] + ':' + str(ln), '|', concl)
