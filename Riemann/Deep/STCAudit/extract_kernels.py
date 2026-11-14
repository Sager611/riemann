
import re, json
L = open('/Users/adriansager/playground/riemann/.scratch/r462/txt/23866074.txt').read().split('\n')
out = {}
i = 0
while i < len(L):
    m = re.match(r'^line_([A-Z]\d+)/sol.json', L[i].strip())
    if m:
        tag = m.group(1)
        # collect until a line that is exactly '}'
        buf = []
        j = i+1
        while j < len(L) and L[j].strip() != '}':
            buf.append(L[j]); j += 1
        buf.append('}')
        txt = '\n'.join(buf)
        # drop page header/footer lines (non-json)
        clean = [ln for ln in txt.split('\n') if not re.search(r'EXPLICIT AND ASYMPTOTIC|CASA|^\s*\d+\s*$', ln)]
        txt = '\n'.join(clean)
        try:
            d = json.loads(txt)
        except Exception as e:
            print(tag, 'PARSE FAIL', e); print(txt[:500]); i = j; continue
        out[tag] = {'d': d['d'], 'B': d['B'], 'c': d['c'], 'm': d['m'], 'C1': d.get('C1'), 'Cz_est': d.get('Cz_est'), 'z_est': d.get('z_est'), 'A3_est': d.get('A3_est'), 'Lstar': d.get('Lstar')}
        i = j
    i += 1
json.dump(out, open('/Users/adriansager/playground/riemann/Riemann/Deep/STCAudit/kernels.json','w'), indent=1)
import math
for k,v in out.items():
    cR1 = v['c'][v['B'].index(['R',1])]
    print(k, 'd', v['d'], 'm', v['m'], 'n', len(v['c']), 'cR1/2pi', cR1/(2*math.pi), 'C1', v['C1'], 'Cz_est', v['Cz_est'], 'z_est', v['z_est'], 'A3_est', v['A3_est'])
