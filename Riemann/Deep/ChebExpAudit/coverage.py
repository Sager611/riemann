
import json, math
from fractions import Fraction as Fr
c = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ChebExp/chains.json'))
tot = 0
for ch in c:
    ivs = ch['intervals']; tot += len(ivs)
    L0 = ch['Lam0']; prev = L0; issues = []
    for k, iv in enumerate(ivs):
        a = iv['Lam_a']
        if Fr(a) > Fr(prev): issues.append(('GAP', k, a, prev, float(Fr(a)-Fr(prev))))
        if iv['kind'] == 'int':
            if not (Fr(iv['Lam_b']) > Fr(a)): issues.append(('EMPTY', k))
            prev = iv['Lam_b']
        else:
            if k != len(ivs)-1: issues.append(('TAIL_NOT_LAST', k))
            prev = math.inf
    kinds = [iv['kind'] for iv in ivs]
    overl = sum(1 for k in range(1,len(ivs)) if ivs[k]['kind']=='int' and Fr(ivs[k]['Lam_a']) < Fr(ivs[k-1]['Lam_b']))
    print(ch['name'], 'c=',ch['c'],'L0=',repr(L0),'first a=',repr(ivs[0]['Lam_a']), 'n=',len(ivs), 'last=',kinds[-1], 'tailstart=',ivs[-1]['Lam_a'], 'strict overlaps=',overl, 'issues=',issues[:5], 'h range', round(min(i['h'] for i in ivs),3), round(max(i['h'] for i in ivs),3), 's0 range', round(min(i['s0'] for i in ivs),4), round(max(i['s0'] for i in ivs),4))
print('total kernels', tot)
print('log13=%.12f log258=%.12f' % (math.log(13), math.log(258)))
