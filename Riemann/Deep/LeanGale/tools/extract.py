#!/usr/bin/env python3
"""Print Lean declarations (with their docstrings) by name: extract.py FILE NAME [NAME ...]"""
import re, sys
path, names = sys.argv[1], sys.argv[2:]
lines = open(path).read().split('\n')
hdr = re.compile(r'^(?:@\[[^\]]*\]\s*)?(?:private |protected )?(?:noncomputable )?(?:def|theorem|lemma|abbrev|structure|inductive|instance) ([^\s:({\[]+)')
for name in names:
    for i, l in enumerate(lines):
        m = hdr.match(l)
        if m and m.group(1) == name:
            s = i
            if s > 0 and lines[s - 1].rstrip().endswith('-/'):
                s -= 1
                while not lines[s].lstrip().startswith('/--'):
                    s -= 1
            e = i + 1
            while e < len(lines) and (lines[e].startswith(' ') or lines[e].startswith('|')) :
                e += 1
            print('\n'.join(lines[s:e]))
            print()
            break
    else:
        print('<<NOT FOUND %s in %s>>' % (name, path))
