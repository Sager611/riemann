# Extracts Tables 2-6 of arXiv 2405.12545v2 (final.tex l.1135-1569) into tables.json (TeX deleted after use; refetch: arxiv.org/e-print/2405.12545)

import re, json
L = open('/Users/adriansager/playground/riemann/.scratch/bellotti2/tex/final.tex').read().split('\n')
def rows(a, b):
    out = []
    for ln in L[a-1:b]:
        if '&' not in ln or 'T_0' in ln: continue
        cells = [c.strip().strip('$').strip() for c in ln.replace('\\\\','').split('&')]
        out.append(cells)
    return out
tabs = {'CB': rows(1135,1220), 'ab': rows(1222,1307), 'uxvw': rows(1309,1394), 'd': rows(1397,1482), 'c': rows(1484,1569)}
for k,v in tabs.items():
    print(k, len(v), v[0], v[-1])
json.dump(tabs, open('/Users/adriansager/playground/riemann/.scratch/bellotti2/tables.json','w'))

