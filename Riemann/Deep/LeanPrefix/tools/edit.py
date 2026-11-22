#!/usr/bin/env python3
"""Exact-text edits of files: python3 edit.py spec.json, spec = [{"file", "old", "new"}, ...].
Each old must occur exactly once in its file (checked before writing)."""
import json, sys
spec = json.load(open(sys.argv[1]))
files = {}
for e in spec:
    f = e['file']
    if f not in files:
        files[f] = open(f).read()
    s = files[f]
    c = s.count(e['old'])
    if c != 1:
        sys.exit('edit.py: %s: old text occurs %d times: %r' % (f, c, e['old'][:80]))
    files[f] = s.replace(e['old'], e['new'])
for f, s in files.items():
    open(f, 'w').write(s)
print('edit.py: ok', ', '.join(files))

