import sys, json, ast; sys.path.insert(0,'.')
from exact4 import check4
import crit3 as C
# consistency with Lemma 1 (J=3): crit3 best point at d24 k3.32: alpha=1, gam, a3, mu -> a=(-1, gam+a3 mu^2, -2 a3 mu, a3), phi12 = ? 
# Lemma-1 form -|f|^2+gam|f'|^2+a3|f'''+mu f'|^2 = -|f|^2+(gam+a3 mu^2)|f'|^2 + a3|f'''|^2 + 2 a3 mu int f''' f'
# and int f''' f' = [f'' f'] - |f''|^2  ->  a2 = -2 a3 mu, phi12 = 2 a3 mu
a3=0.00132656; mu=4.85647; gam=0.98659291
print('J3 via exact4:', check4([-1, gam+a3*mu*mu, -2*a3*mu, a3], 2*a3*mu, 0.0, ng=200, KM=4000))
print('J3 via crit3 valid():', C.valid(gam, a3, mu))
for fn in ('logs/scanJ_r2_d24.json','logs/scanJ_r2_d64.json'):
    for line in json.load(open(fn)):
        a=line['a4']; ph={ast.literal_eval(k):v for k,v in line['phi4'].items()}
        print(fn[-9:], line['kappa'], 'h4=%.6f'%line['h4'], check4(a, ph[(1,2)], ph[(3,2)], ng=300, KM=20000), flush=True)
