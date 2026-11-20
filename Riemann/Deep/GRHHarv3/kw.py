import re, sys
# usage: python3 kw.py FILE kw1 kw2 ... : first hit of each keyword with 150/180 chars of context
f = sys.argv[1]
t = open(f, errors='replace').read()
for k in sys.argv[2:]:
    ms = list(re.finditer(k, t, re.I))
    if not ms:
        print(k, '| none'); continue
    m = ms[0]
    print(k, '| n=%d | l.%d |' % (len(ms), t.count('\n', 0, m.start()) + 1),
          t[max(0, m.start() - 150):m.start() + 180].replace('\n', ' '))
